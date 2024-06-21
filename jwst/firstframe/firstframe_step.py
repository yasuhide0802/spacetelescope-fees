from ..stpipe import Step
from . import firstframe_sub
from jwst.lib.basic_utils import use_datamodel


__all__ = ["FirstFrameStep"]


class FirstFrameStep(Step):
    """
    FirstFrameStep: This is a MIRI specific task.  If the number of groups
    is greater than 3, the DO_NOT_USE group data quality flag is added to
    first group.
    """

    class_alias = "firstframe"

    def process(self, input):

        # Open the input data model
        input_model = use_datamodel(input)

        # check the data is MIRI data
        detector = input_model.meta.instrument.detector.upper()
        if detector[:3] == 'MIR':
            # Do the firstframe correction subtraction
            result = firstframe_sub.do_correction(input_model)
        else:
            self.log.warning('First Frame Correction is only for MIRI data')
            self.log.warning('First frame step will be skipped')
            result = input_model.copy()
            result.meta.cal_step.firstframe = 'SKIPPED'

        input_model.close()
        return result
