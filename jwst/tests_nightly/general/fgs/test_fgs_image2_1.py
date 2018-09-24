import pytest
from jwst.pipeline.calwebb_image2 import Image2Pipeline

from jwst.tests.base_test import BaseJWSTTest


@pytest.mark.bigdata
class TestImage2Pipeline(BaseJWSTTest):
    input_loc = 'fgs'
    ref_loc = ['test_image2pipeline', 'truth']

    def test_fgs_image2pipeline1(self, _jail):
        """

        Regression test of calwebb_image2 pipeline performed on FGS imaging mode data.

        """
        input_file = self.get_data('test_image2pipeline',
                                   'jw86500007001_02101_00001_GUIDER2_rate.fits')
        output_file = 'jw86500007001_02101_00001_GUIDER2_cal.fits'

        Image2Pipeline.call(input_file,output_file = output_file)

        outputs = [(output_file, 'jw86500007001_02101_00001_GUIDER2_cal_ref.fits',
                    ['primary','sci','err','dq','area'])
                  ]
        self.compare_outputs(outputs)
