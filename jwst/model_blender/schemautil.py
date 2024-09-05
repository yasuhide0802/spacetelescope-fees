from stdatamodels import schema as dm_schema

def parse_schema(schema):
    """
    Returns
    -------
    attr_to_column: dict
        one-to-one mapping of metadata attributes
        (full dotted-paths) to table columns. For example
        {'meta.observation.time': 'TIME-OBS'} for storing
        the 'meta.observation.time' attributes in a 'TIME-OBS'
        column.

    attr_to_rule: dict
        mapping of metadata attribute to blend rule
        type. For example {'meta.observation.time': 'mintime'}
        will combine all 'meta.observation.time' attributes
        using the 'mintime' rule.
    """
    def callback(subschema, path, combiner, ctx, recurse):
        if len(path) <= 1:
            return  # ignore top-level (non-meta) attributes
        if path[0] != 'meta':
            return  # ignore non-meta attributes
        if 'items' in path:
            return  # ignore attributes in arrays
        if subschema.get('type') == 'array':
            return  # ignore ListNodes
        if subschema.get('properties'):
            return  # ignore ObjectNodes

        # strip trailing path if there's a combiner
        for combiner in ['anyOf', 'oneOf']:
            if combiner in path:
                path = path[:path.index(combiner)]
                break

        # construct the metadata attribute path
        attr = '.'.join(path)

        # if 'blend_rule' is defined, make a 'blend'
        if 'blend_rule' in subschema:
            ctx['blends'][attr] = subschema['blend_rule']

        # if 'blend_table' is defined (and truthy), add a column
        if subschema.get('blend_table'):
            ctx['columns'][attr] = subschema.get('fits_keyword', attr)

    ctx = {
        'columns': {},
        'blends': {},
    }
    dm_schema.walk_schema(schema, callback, ctx)
    return ctx['columns'], ctx['blends']
