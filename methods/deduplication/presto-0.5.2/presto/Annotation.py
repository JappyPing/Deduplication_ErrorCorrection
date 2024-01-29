"""
Annotation functions
"""
# Info
__author__ = 'Jason Anthony Vander Heiden'
from presto import __version__, __date__

# Imports
from collections import OrderedDict
from presto.Defaults import default_delimiter, default_coord_type


def parseAnnotation(record, fields=None, delimiter=default_delimiter):
    """
    Extracts annotations from a FASTA/FASTQ sequence description

    Arguments:
      record : Description string to extract annotations from
      fields : List of fields to subset the return dictionary to;
               if None return all fields
      delimiter : a tuple of delimiters for (fields, values, value lists)

    Returns:
      OrderedDict : An OrderedDict of field/value pairs
    """
    annotation = record.split(delimiter[0])
    field_dict = OrderedDict([('ID', annotation.pop(0))])
    for ann in annotation:
        vals = ann.split(delimiter[1])
        field_dict[vals[0].upper()] = vals[1]

    # Subset return dictionary to requested fields
    if fields is not None:
        if not isinstance(fields, list):  fields = [fields]
        for f in set(field_dict).difference(fields):  del field_dict[f]

    return field_dict


def flattenAnnotation(ann_dict, delimiter=default_delimiter):
    """
    Converts annotations from a dictionary to a FASTA/FASTQ sequence description

    Arguments:
      ann_dict : Dictionary of field/value pairs
      delimiter : Tuple of delimiters for (fields, values, value lists)

    Returns:
      str : Formatted sequence description string
    """
    annotation = ann_dict.get('ID', 'NONE')
    for k, v in ann_dict.items():
        # Skip ID field
        if k.upper() == 'ID':
            continue

        if isinstance(v, list):
            v = delimiter[2].join([str(x) for x in v])
        annotation += '%s%s%s%s' % (delimiter[0], k.upper(), delimiter[1], v)

    return annotation


def mergeAnnotation(ann_dict_1, ann_dict_2, prepend=False,
                    delimiter=default_delimiter):
    """
    Merges non-ID field annotations from one field dictionary into another

    Arguments:
      ann_dict_1 : Dictionary of field/value pairs to append to
      ann_dict_2 : Dictionary of field/value pairs to merge with ann_dict_2
      prepend : If True then add ann_dict_2 values to the front of any ann_dict_1
                values that are already present, rather than the default behavior
                of appending ann_dict_2 values.
      delimiter : Tuple of delimiters for (fields, values, value lists)

    Returns:
      OrderedDict : Modified ann_dict_1 dictonary of field/value pairs
    """
    # Define merge order
    if prepend:
        def _merge(x, y):  return '%s%s%s' % (y, delimiter[2], x)
    else:
        def _merge(x, y):  return '%s%s%s' % (x, delimiter[2], y)

    merged_dict = ann_dict_1.copy()
    for k, v in ann_dict_2.items():
        # Skip ID field
        if k.upper() == 'ID':
            continue

        if k in merged_dict:
            if isinstance(merged_dict[k], list):
                merged_dict[k] = delimiter[2].join([str(x) for x in merged_dict[k]])
            if isinstance(v, list):
                v = delimiter[2].join([str(x) for x in v])
            merged_dict[k] = _merge(merged_dict[k], v)
        else:
            merged_dict[k.upper()] = v

    return merged_dict


def renameAnnotation(ann_dict, old_field, new_field, delimiter=default_delimiter):
    """
    Renames an annotation and merges annotations if the new name already exists

    Arguments:
      ann_dict : Dictionary of field/value pairs
      old_field : Old field name
      new_field : New field name
      delimiter : Tuple of delimiters for (fields, values, value lists)

    Returns:
      OrderedDict : Modified fields dictonary
    """
    if new_field in ann_dict:
        rename_dict = ann_dict.copy()
        del rename_dict[old_field]
        rename_dict = mergeAnnotation(rename_dict, {new_field:ann_dict[old_field]},
                                      delimiter=delimiter)
    else:
        rename_dict = OrderedDict([(new_field, v) if k == old_field else (k, v) \
                                   for k, v in ann_dict.items()])

    return rename_dict


# TODO:  this converted min/max/sum collapse to strings instead of floats (for rounding purposes). which is odd.
def collapseAnnotation(ann_dict, action, fields=None, delimiter=default_delimiter):
    """
    Collapses multiple annotations into new single annotations for each field

    Arguments:
      ann_dict : Dictionary of field/value pairs
      action : Collapse action to take;
               one of {min, max, sum, first, last, set, cat}
      fields : Subset of ann_dict to _collapse;
               if None _collapse all but the ID field
      delimiter : Tuple of delimiters for (fields, values, value lists)

    Returns:
      OrderedDict : Modified field dictionary
    """
    # Define _collapse action
    if action == 'set':
        def _collapse(value):  return sorted(set(value))
    elif action == 'first':
        def _collapse(value):  return value[0]
    elif action == 'last':
        def _collapse(value):  return value[-1]
    elif action == 'min':
        def _collapse(value):  return '%.12g' % min([float(x or 0) for x in value])
    elif action == 'max':
        def _collapse(value):  return '%.12g' % max([float(x or 0) for x in value])
    elif action == 'sum':
        def _collapse(value):  return '%.12g' % sum([float(x or 0) for x in value])
    elif action == 'cat':
        def _collapse(value):  return ''.join([str(x) for x in value])
    else:
        def _collapse(value):  return value

    # Collapse fields
    collapse_dict = ann_dict.copy()
    for k, v in collapse_dict.items():
        if k.upper() == 'ID':
            continue
        if fields is None or k in fields:
            # Convert field to list
            if not isinstance(v, list) and isinstance(v, str):
                v = v.split(delimiter[2])
            elif not isinstance(v, list):
                v = [v]
            # Perform _collapse and reassign field
            collapse_dict[k] = _collapse(v)

    return collapse_dict


def getAnnotationValues(seq_iter, field, unique=False, delimiter=default_delimiter):
    """
    Gets the set of unique annotation values in a sequence set

    Arguments:
      seq_iter : Iterator or list of SeqRecord objects
      field : Annotation field to retrieve values for
      unique : If True return a list of only the unique values;
               if False return a list of all values
      delimiter : Tuple of delimiters for (fields, values, value lists)

    Returns:
      list : List of values for the field
    """
    # Parse annotations from seq_list records
    ann_iter = (parseAnnotation(s.description, delimiter=delimiter) for s in seq_iter)
    values = [a[field] for a in ann_iter]

    return list(set(values)) if unique else values


def annotationConsensus(seq_iter, field, delimiter=default_delimiter):
    """
    Calculate a consensus annotation for a set of sequences

    Arguments:
      seq_iter : an iterator or list of SeqRecord objects
      field : the annotation field to take a consensus of
      delimiter : a tuple of delimiters for (annotations, field/values, value lists)

    Returns:
      dict : Dictionary with keys
             `set` containing a list of unique annotation values,
             `count` containing annotation counts,
             `cons` containing the consensus annotation,
             `freq` containing the majority annotation frequency
    """
    # Define return dictionary
    cons_dict = {'set':None, 'count':None, 'cons':None, 'freq':None}

    # Parse annotations from seq_list records
    val_list = getAnnotationValues(seq_iter, field, delimiter=delimiter)

    # Define annotation set and counts
    cons_dict['set'] = sorted(set(val_list))
    cons_dict['count'] = [val_list.count(v) for v in cons_dict['set']]

    # Define consensus annotation
    i = cons_dict['count'].index(max(cons_dict['count']))
    cons_dict['cons'] = cons_dict['set'][i]
    cons_dict['freq'] = float(cons_dict['count'][i]) / len(val_list)

    return cons_dict


def getCoordKey(header, coord_type=default_coord_type, delimiter=default_delimiter):
    """
    Return the coordinate identifier for a sequence description

    Arguments:
      header : Sequence header string
      coord_type : Sequence header format;
                   one of ['illumina', 'solexa', 'sra', '454', 'presto'];
                   if unrecognized type or None return sequence ID.
      delimiter : Tuple of delimiters for (fields, values, value lists)

    Returns:
      str : Coordinate identifier as a string
    """
    #header = seq.id
    if coord_type in ('illumina', 'solexa'):
        return header.split()[0].split('#')[0]
    elif coord_type == '454':
        return header.split()[0]
    elif coord_type == 'sra':
        return '.'.join(header.split()[0].split('.')[:2])
    elif coord_type == 'presto':
        return parseAnnotation(header, delimiter=delimiter)['ID']
    else:
        return header
