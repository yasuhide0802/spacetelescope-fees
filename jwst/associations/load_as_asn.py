"""Treat various objects as Associations"""

from os.path import splitext

from ..associations import (
    Association,
    AssociationRegistry,
    libpath,
    load_asn
)
from ..associations.asn_from_list import asn_from_list
from ..associations.lib.rules_level2_base import DMSLevel2bBase

__all__ = [
    'LoadAsAssociation',
    'LoadAsLevel2Asn',
]


DEFAULT_NAME = 'cast_to_association'
DEFAULT_ASN_META = {
    'program': DEFAULT_NAME,
    'target': DEFAULT_NAME,
    'asn_pool': DEFAULT_NAME
}


class LoadAsAssociation(dict):
    """Read in or create an association

    Parameters
    ----------
    asn: dict or Association
        An already existing association

    Notes
    -----
    This class is normally not instantiated.
    the `load` method should be used as the factory
    method to read an association or create one from
    a string or `Datamodel` object, or a list of such
    objects.
    """

    @classmethod
    def load(cls, obj,
             meta=DEFAULT_ASN_META,
             registry=AssociationRegistry,
             rule=Association,
             product_name_func=None):
        """ Load object and return an association of it

        Parameters
        ----------
        obj: Association, str, Datamodel, [str[,...]], [Datamodel[,...]]
            The obj to return as an association

        registry: AssociationRegistry
            The registry to use to load an association file with

        rule: Association
            The rule to use if an association needs to be created

        product_name_func: func
            A function, when given the argument of `obj`, or
            if `obj` is a list, each item in `obj`, returns
            a string that will be used as the product name in
            the association.

        Attributes
        ----------
        Along with the attributes belonging to an association, the
        following are added:

        filename: str
            The name of the association file, if such a file
            where passed in. Otherwise a default value is given.

        Returns
        -------
        association: Association
            An association created using given obj
        """
        try:
            with open(obj) as fp:
                pure_asn = load_asn(fp, registry=registry)
        except Exception:
            if not isinstance(obj, list):
                obj = [obj]
            asn = asn_from_list(
                obj,
                rule=rule,
                meta=DEFAULT_ASN_META,
                product_name_func=product_name_func
            )
            asn.filename = DEFAULT_NAME
        else:
            asn = cls(pure_asn)
            asn.filename = obj

        return asn


class LoadAsLevel2Asn(LoadAsAssociation):
    """Read in or create a Level2 association
    """

    @classmethod
    def load(cls, obj):
        """ Open object and return a Level2 association of it

        Parameters
        ----------
        obj: Association, str, Datamodel, [str[,...]], [Datamodel[,...]]
            The obj to return as an association

        Attributes
        ----------
        Along with the attributes belonging to a Level2 association, the
        following are added:

        filename: str
            The name of the association file, if such a file
            where passed in. Otherwise a default value is given.

        Returns
        -------
        association: DMSLevel2bBase
            An association created using given obj
        """
        asn = super(LoadAsLevel2Asn, cls).load(
            obj,
            registry=AssociationRegistry(
                definition_files=[libpath('rules_level2b.py')],
                include_default=False
            ),
            rule=DMSLevel2bBase,
            product_name_func=cls.model_product_name
        )
        return asn

    @staticmethod
    def model_product_name(model):
        """Product a model product name based on the model.

        Parameters
        ----------
        model: DataModel
            The model to get the name from
        """
        product_name, ext = splitext(model.meta.filename)
        return product_name
