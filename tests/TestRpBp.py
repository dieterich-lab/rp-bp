###
#   This class contains unit and (eventually) integration tests for rpbp.
#   The structure of this file is taken from an old blog post here:
#       http://blog.jameskyle.org/2010/10/nose-unit-testing-quick-start/   
###

import os

import appdirs

class TestRpBp(object):
    
    @classmethod
    def prepare_data(self):
        # first, download the data to <user_data_dir>/rpbp/test-data
        loc = os.path.join(appdirs.user_data_dir("rpbp"), "test-data")
        

        pass

    @classmethod
    def setup_class(klass):
        """This method is run once for each class before any tests are run"""
        TestRpBp.prepare_data()
        

    @classmethod
    def teardown_class(klass):
        """This method is run once for each class _after_ all tests are run"""
        pass

    def setUp(self):
        """This method is run once before _each_ test method is executed"""
        pass

    def teardown(self):
        """This method is run once after _each_ test method is executed"""
        pass


