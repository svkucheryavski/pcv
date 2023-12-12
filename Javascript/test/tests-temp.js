/*****************************************************************
 *  Tests for PCV methods                                        *
 *****************************************************************/

// import dependencies
import {default as chai} from 'chai';
import {default as chaiAlmost} from 'chai-almost';
import { matrix, vector, index, c, Index } from 'mdatools/arrays';
import { pcafit, pcapredict, pcrfit, pcrpredict, plsfit, plspredict, splitregdata } from 'mdatools/models';

// import methods to test
import { cv2obs, crossval, pcvpca, pcvpcr, pcvpls } from '../src/index.js';

// set up test settings
const expect = chai.expect;
chai.use(chaiAlmost(0.0001));

describe('Manual tests for PCV methods', function () {

   // common dataset for PCR and PLS tests
   const data = matrix([
      32, 150, 41, 28000, 119,
      35, 160, 48, 31000, 129,
      36, 166, 47, 28000, 112,
      37, 166, 49, 14000, 123,
      42, 175, 67, 38000, 105,
      43, 180, 80, 30000, 129,
      43, 181, 75, 31000, 105,
      44, 180, 81, 42000, 113,
   ], 5, 8).t();


   it ('tests for method "pcvpca"', function () {

      // comparing outcomes with R
      const X1 = matrix([1, 3, 17, 19, 10, 14, 7, 13, 9, 11, 2, 15, 8, 6, 4, 12,
         22, 18, 4, 14, 13, 12, 1, 15], 8, 3);


      // second case - 2 components + ven with 3 segments
      const A12 = 2;
      const cv12 = {type: "ven", nseg: 3};
      const m12 = pcafit(X1, A12);
      const X12pv = pcvpca(X1, m12, A12, cv12);
   });

});
