# coding=utf-8

import unittest
import pySOS


class TestOSOAAClasses(unittest.TestCase):

        def testAER(self):
            aer = pySOS.AER()
            self.assertEqual(aer.waref, 0.55)
            self.assertEqual(aer.aotref, 0.1)
            self.assertIsNone(aer.tronca)
            self.assertEqual(aer.model, 2)
            self.assertEqual(aer.sf.model, 3)
            self.assertEqual(aer.sf.rh, 98)

            aer.SetModel(model=0, sdtype=1)
            self.assertIsNone(aer.wmo)
            self.assertIsNone(aer.sf)
            self.assertIsNone(aer.lnd)
            self.assertIsNone(aer.external)
            self.assertEqual(aer.model, 0)
            self.assertEqual(aer.mm.sdtype, 1)
            self.assertListEqual(list(vars(aer.mm).keys()), ["sdtype",
                                                             "mrwa", "miwa",
                                                             "mrwaref", "miwaref",
                                                             "sdradius", "sdvar"])

            aer.SetModel(model=0, sdtype=2)
            self.assertIsNone(aer.wmo)
            self.assertIsNone(aer.sf)
            self.assertIsNone(aer.lnd)
            self.assertIsNone(aer.external)
            self.assertEqual(aer.model, 0)
            self.assertEqual(aer.mm.sdtype, 2)
            self.assertListEqual(list(vars(aer.mm).keys()), ["sdtype", "mrwa",
                                                             "miwa", "mrwaref",
                                                             "miwaref", "slope",
                                                             "rmin", "rmax"])

            aer.SetModel(model=1, wmotype=2)
            self.assertIsNone(aer.mm)
            self.assertIsNone(aer.sf)
            self.assertIsNone(aer.lnd)
            self.assertIsNone(aer.external)
            self.assertEqual(aer.model, 1)
            self.assertEqual(aer.wmo.model, 2)

            aer.SetModel(model=1, wmotype=4, dl=0.4, ws=0.3, oc=0.2, so=0.1)
            self.assertIsNone(aer.mm)
            self.assertIsNone(aer.sf)
            self.assertIsNone(aer.lnd)
            self.assertIsNone(aer.external)
            self.assertEqual(aer.model, 1)
            self.assertEqual(aer.wmo.model, 4)
            self.assertEqual(aer.wmo.dl, 0.4)
            self.assertEqual(aer.wmo.ws, 0.3)
            self.assertEqual(aer.wmo.oc, 0.2)
            self.assertEqual(aer.wmo.so, 0.1)

            aer.SetModel(model=2, sfmodel=2, rh=50)
            self.assertEqual(aer.model, 2)
            self.assertEqual(aer.sf.model, 2)
            self.assertEqual(aer.sf.rh, 50)
            self.assertIsNone(aer.mm)
            self.assertIsNone(aer.wmo)
            self.assertIsNone(aer.lnd)
            self.assertIsNone(aer.external)

            aer.SetModel(model=3, vcdef=1)
            self.assertEqual(aer.model, 3)
            self.assertEqual(aer.lnb.vcdef, 1)
            self.assertIsNone(aer.mm)
            self.assertIsNone(aer.wmo)
            self.assertIsNone(aer.sf)
            self.assertIsNone(aer.external)
            self.assertListEqual(list(vars(aer.lnb).keys()), ["vcdef",
                                                              "coarsevc", "finevc",
                                                              "cmrwa", "cmiwa",
                                                              "cmrwaref",
                                                              "cmiwaref",
                                                              "csdradius",
                                                              "csdvar",
                                                              "fmrwa", "fmiwa",
                                                              "fmrwaref",
                                                              "fmiwaref",
                                                              "fsdradius",
                                                              "fsdvar"])

            aer.SetModel(model=3, vcdef=2)
            self.assertEqual(aer.model, 3)
            self.assertEqual(aer.lnb.vcdef, 2)
            self.assertIsNone(aer.mm)
            self.assertIsNone(aer.wmo)
            self.assertIsNone(aer.sf)
            self.assertIsNone(aer.external)
            self.assertListEqual(list(vars(aer.lnb).keys()), ["vcdef",
                                                              "raot",
                                                              "cmrwa", "cmiwa",
                                                              "cmrwaref",
                                                              "cmiwaref",
                                                              "csdradius",
                                                              "csdvar",
                                                              "fmrwa", "fmiwa",
                                                              "fmrwaref",
                                                              "fmiwaref",
                                                              "fsdradius",
                                                              "fsdvar"])

            aer.SetModel(model=4, extdata="test.txt")
            self.assertEqual(aer.model, 4)
            self.assertEqual(aer.extdata, "test.txt")
            self.assertIsNone(aer.mm)
            self.assertIsNone(aer.wmo)
            self.assertIsNone(aer.sf)
            self.assertIsNone(aer.lnd)
