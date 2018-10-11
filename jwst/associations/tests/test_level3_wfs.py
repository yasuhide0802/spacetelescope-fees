"""test_level3_dithers: Test of WFS rules."""

from . import helpers

from .. import generate

# Generate Level3 assocations
rules = helpers.registry_level3_only()
pool = helpers.combine_pools(
    helpers.t_path('data/pool_004_wfs.csv')
)
level3_asns = generate(pool, rules)


class TestLevel3WFS(helpers.BasePoolRule):

    pools = [
        helpers.PoolParams(
            path=helpers.t_path('data/pool_004_wfs.csv'),
            n_asns=35,
            n_orphaned=0
        ),
    ]

    valid_rules = [
        'Asn_WFSCMB',
    ]


def test_wfs_duplicate_product_names():
    """Test for duplicate product names"""
    global level3_asns

    name_list = [
        product['name']
        for asn in level3_asns
        for product in asn['products']
    ]
    assert len(name_list)
    name_set = set(name_list)
    assert len(name_set) == len(name_list)
