"""Read information from config dictionaries

Functions:
  - read_param: Read a parameter from a dictionary
"""

import logging
logger = logging.getLogger(__name__)

def read_param(yaml_data, node, default):
    """Read a parameter from a dictionary

    Args:
        yaml_data: dictionary containing data to read
        node: Path through dictionaries to read. For example,
            node = "key1.key2.key3" will try to return
            yaml_data["key1"]["key2"]["key3"]
        default: Value to return if the given node cannot be found.
            if default=='required', then an exception will be thrown
            if the node cannot be found
    """
    data = yaml_data
    xpath = node.split('.')
    try:
        for path in xpath:
            data = data[path]
    except Exception as exc:
        if default == 'required':
            err = "FATAL: Configuration parameter " + \
                  """{0}""".format(node) + \
                  " required but not provided in config file!"
            logger.debug(err)
            raise exc
        else:
            data = default
    return data
