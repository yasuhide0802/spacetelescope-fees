import logging

from collections import deque
import numpy as np

from .association import (
    ProcessList,
    make_timestamp
)

# Configure logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def generate(pool, rules, version_id=None):
    """Generate associations in the pool according to the rules.

    Parameters
    ----------
    pool: AssociationPool
        The pool to generate from.

    rules: Associations
        The associaton rule set.

    version_id: None, True, or str
        The string to use to tag associations and products.
        If None, no tagging occurs.
        If True, use a timestamp
        If a string, the string.

    Returns
    -------
    associations: [association[,...]]
        List of associations

    Notes
    -----
    Refer to the :ref:`Association Generator <association-generator>`
    documentation for a full description.
    """
    associations = []
    in_an_asn = np.zeros((len(pool),), dtype=bool)
    if type(version_id) is bool:
        version_id = make_timestamp()
    process_list = ProcessQueue([
        ProcessList(
            items=pool,
            rules=[rule for _, rule in rules.items()]
        )
    ])

    for process_idx, process_item in enumerate(process_list):
        for item in process_item.items:

            # Determine against what the item should be compared
            # against.
            use_rules = rules
            use_associations = associations
            if process_item.work_over == process_item.EXISTING:
                use_rules = None
            if process_item.work_over == process_item.RULES:
                use_associations = []
            existing_asns, new_asns, to_process = generate_from_item(
                item,
                version_id,
                use_associations,
                use_rules,
                process_item.rules
            )
            associations.extend(new_asns)

            # If working on a process list EXISTING
            # remove any new `to_process` that is
            # also EXISTING. Prevent infinite loops.
            if process_item.work_over == ProcessList.EXISTING:
                to_process = [
                    to_process_item
                    for to_process_item in to_process
                    if to_process_item.work_over != ProcessList.EXISTING
                ]
            process_list.extend(to_process)

            if len(existing_asns) +\
               len(new_asns) > 0:
                in_an_asn[item.index] = True

    # Finalize found associations
    finalized_asns = rules.finalize(associations)

    return finalized_asns


def generate_from_item(
        item,
        version_id,
        associations,
        rules,
        allowed_rules):
    """Either match or generate a new assocation

    Parameters
    ----------
    item: dict
        The item to match to existing associations
        or generate new associations from

    version_id: str or None
        Version id to use with association creation.
        If None, no versioning is used.

    associations: [association, ...]
        List of already existing associations.
        If the item matches any of these, it will be added
        to them.

    rules: AssociationRegistry or None
        List of rules to create new associations

    allowed_rules: [rule, ...]
        Only compare existing associations and rules if the
        rule is in this list. If none,
        all existing associations and rules will be checked.

    Returns
    -------
    (associations, process_list): 3-tuple where
        existing_asns: [association,...]
            List of existing associations item belongs to.
            Empty if none match
        new_asns: [association,...]
            List of new associations item creates. Empty if none match
        process_list: [ProcessList, ...]
            List of process events.
    """

    # Check membership in existing associations.
    associations = [
        asn
        for asn in associations
        if type(asn) in allowed_rules
    ]
    existing_asns, process_list = match_item(item, associations)

    # Now see if this item will create new associatons.
    # By default, a item will not be allowed to create
    # an association based on rules of existing associations.
    to_process = []
    new_asns = []
    if rules is not None:
        ignore_asns = set([type(asn) for asn in existing_asns])
        new_asns, to_process = rules.match(
            item,
            version_id=version_id,
            allow=allowed_rules,
            ignore=ignore_asns,
        )
    process_list.extend(to_process)

    return existing_asns, new_asns, process_list


def match_item(item, associations):
    """Match item to a list of associations

    Parameters
    ----------
    item: dict
        The item to match to the associations.

    associations: [association, ...]
        List of already existing associations.
        If the item matches any of these, it will be added
        to them.

    Returns
    -------
    (associations, process_list): 2-tuple where
        associations: [association,...]
            List of associations item belongs to. Empty if none match
        process_list: [ProcessList, ...]
            List of process events.
    """
    item_associations = []
    process_list = []
    for asn in associations:
        if asn in item_associations:
            continue
        matches, reprocess = asn.add(item)
        process_list.extend(reprocess)
        if matches:
            item_associations.append(asn)
    return item_associations, process_list


class ProcessQueue(deque):
    """Make a deque iterable and mutable"""
    def __iter__(self):
        while True:
            try:
                yield self.popleft()
            except:
                break
