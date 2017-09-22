from __future__ import (division, print_function, unicode_literals,
    absolute_import)

from ..stpipe import Step, cmdline
from .. import datamodels
from . import resample_spec
from ..exp_to_source import multislit_to_container


class ResampleSpecStep(Step):
    """
    ResampleSpecStep: Resample input data onto a regular grid using the
    drizzle algorithm.

    Parameters
    -----------
    input : DataModel, Association
    """

    spec = """
        single = boolean(default=False)
        wht_type = option('exptime', 'error', None, default='exptime')
        pixfrac = float(default=1.0)
        kernel = string(default='square')
        fillval = string(default='INDEF')
        good_bits = integer(default=4)
    """
    reference_file_types = ['drizpars']

    def process(self, input):

        input_models = datamodels.open(input)

        # Put single input into a ModelContainer
        if not isinstance(input_models, datamodels.ModelContainer):
            s = datamodels.ModelContainer()
            s.append(input_models)
            input_models = s

        self.driz_filename = self.get_reference_file(input_models[0], 'drizpars')

        # Multislits get converted to a ModelContainer per slit
        if all([isinstance(i, datamodels.MultiSlitModel) for i in input_models]):
            self.log.info('Converting MultiSlit to ModelContainer')
            container_dict = multislit_to_container(input_models)
            output_product = datamodels.MultiProductModel()
            output_product.update(input_models[0])
            for k, v in container_dict.items():
                input_models = v

                # Set up the resampling object as part of this step
                self.step = resample_spec.ResampleSpecData(input_models,
                    self.driz_filename, single=self.single,
                    wht_type=self.wht_type, pixfrac=self.pixfrac,
                    kernel=self.kernel, fillval=self.fillval,
                    good_bits=self.good_bits)
                # Do the resampling
                self.step.do_drizzle()
                if len(self.step.output_models) == 1:
                    out_slit = self.step.output_models[0]
                    output_product.products.append(out_slit)
                else:
                    out_slit = self.step.output_models
            result = output_product
        else:
            # Set up the resampling object as part of this step
            self.step = resample_spec.ResampleSpecData(input_models,
                self.driz_filename, single=self.single, wht_type=self.wht_type,
                pixfrac=self.pixfrac, kernel=self.kernel,
                fillval=self.fillval, good_bits=self.good_bits)
            # Do the resampling
            self.step.do_drizzle()

            # Return either the single resampled datamodel, or the container
            # of datamodels.
            if len(self.step.output_models) == 1:
                result = self.step.output_models[0]
            else:
                result = self.step.output_models

        result.meta.cal_step.resample = 'COMPLETE'

        return result


if __name__ == '__main__':
    cmdline.step_script(ResampleSpecStep)
