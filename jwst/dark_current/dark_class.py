import numpy as np

class DarkData:
    """
    Class removing data model dependencies.
    """
    def __init__(self, dims=None, dark_model=None):
        if dark_model is not None:
            self.data = dark_model.data
            self.dq = dark_model.dq
            self.err = dark_model.err

            self.exp_nframes = dark_model.meta.exposure.nframes
            self.exp_ngroups = dark_model.meta.exposure.ngroups
            self.exp_groupgap = dark_model.meta.exposure.groupgap

        elif dims is not None:
            self.data = np.zeros(dims, dtype=np.float32)
            self.dq = np.zeros(dims, dtype=np.uint8)
            self.err = np.zeros(dims, dtype=np.float32)

            self.exp_nframes = None
            self.exp_ngroups = None
            self.exp_groupgap = None

        else:
            self.data = None
            self.dq = None
            self.err = None

            self.exp_nframes = None
            self.exp_ngroups = None
            self.exp_groupgap = None

        self.save = False
        self.oname = None


class DarkScienceData:
    """
    Class removing data model dependencies.
    """
    def __init__(self, science_model):
        self.data = science_model.data
        self.dq = science_model.groupdq
        self.pixeldq = science_model.pixeldq
        self.err = science_model.err

        self.instrument_name = science_model.meta.instrument.name
        self.exp_nframes = science_model.meta.exposure.nframes
        self.exp_groupgap = science_model.meta.exposure.groupgap

        self.cal_step = None
