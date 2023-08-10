/********************************************************************
 *  Tests for PCV methods by comparing with reference values from R *
 ********************************************************************/

// import dependencies
import * as fs from 'fs';
import {default as chai} from 'chai';
import {default as chaiAlmost} from 'chai-almost';

import { quantile, max } from 'mdatools/stat';
import { ismatrix, Matrix } from 'mdatools/arrays';
import { pcafit, pcapredict, pcrfit, pcrpredict, plsfit, plspredict, splitregdata } from 'mdatools/models';

// import methods to test
import { pcvpca, pcvpcr, pcvpls } from '../src/index.js';

// set up test settings
const center = true;
const expect = chai.expect;

// directory with reference test values from R and function for checking its existence
const referenceDir = '../.tests/';
function checkTestDir(dir) {
   return fs.existsSync(referenceDir) && fs.existsSync(referenceDir + 'data') && fs.existsSync(referenceDir + dir);
}

// tolerance values to expect from the tests
const tolH = 0.001;  // max difference is 0.1% of the 95th percentile for reference distances
const tolQ = 0.001;  // max difference is 0.1% of the 95th percentile for reference distances
const tolD = 0.010;  // max absolute difference
const tolY = 0.010;  // max absolute difference

chai.use(chaiAlmost());

/**
 * Read reference values from file and return as a matrix.
 * @param {string} path - path to file.
 * @returns {Matrix}
 */
function readOutcome(path) {
   return Matrix.parseCSV(fs.readFileSync(path, 'utf8')).values;
}

/**
 * Runs PCA based test.
 * @param {Matrix} X - matrix with predictors.
 * @param {Matrix} Y - matrix with responses (not needed just a placeholder).
 * @param {ncomp} ncomp - number of components for PV-set creation.
 * @param {boolean} scale - to scale or not to scale calibration set when creating PV-set.
 * @param {JSON} cv - cross-validation resampling parameters.
 * @param {string} fileSuffix - suffix for reference files.
 */
function runPCATest(X, Y, ncomp, scale, cv, fileSuffix) {

   // compute results for global and local CV-scopes
   const m = pcafit(X, ncomp, center, scale);

   const Xpvg = pcvpca(X, m, ncomp, cv, 'global');
   const Xpvl = pcvpca(X, m, ncomp, cv, 'local');

   const rg = pcapredict(m, Xpvg);
   const rl = pcapredict(m, Xpvl);

   // read reference values
   const tHg = readOutcome(referenceDir + 'pcvpca/' + 'Hpvg' + fileSuffix);
   const tQg = readOutcome(referenceDir + 'pcvpca/' + 'Qpvg' + fileSuffix);
   const tHl = readOutcome(referenceDir + 'pcvpca/' + 'Hpvl' + fileSuffix);
   const tQl = readOutcome(referenceDir + 'pcvpca/' + 'Qpvl' + fileSuffix);

   // tests for correct type
   expect(ismatrix(tHg) && ismatrix(rg.H)).to.be.true;
   expect(ismatrix(tQg) && ismatrix(rg.Q)).to.be.true;
   expect(ismatrix(tHl) && ismatrix(rl.H)).to.be.true;
   expect(ismatrix(tQl) && ismatrix(rl.Q)).to.be.true;

   // compute maximum absolute relative difference between distances as % of 95th reference distances percentile
   const hgMaxDiff = max(tHg.subtract(rg.H).apply(x => Math.abs(x), 0).v) / quantile(tHg.v, 0.95)
   const qgMaxDiff = max(tQg.subtract(rg.Q).apply(x => Math.abs(x), 0).v) / quantile(tQg.v, 0.95)
   const hlMaxDiff = max(tHl.subtract(rl.H).apply(x => Math.abs(x), 0).v) / quantile(tHl.v, 0.95)
   const qlMaxDiff = max(tQl.subtract(rl.Q).apply(x => Math.abs(x), 0).v) / quantile(tQl.v, 0.95)

   // test that the difference does not exceed tolerance limits
   expect(hgMaxDiff < tolH).to.be.true;
   expect(qgMaxDiff < tolQ).to.be.true;
   expect(hlMaxDiff < tolH).to.be.true;
   expect(qlMaxDiff < tolQ).to.be.true;
}

/**
 * Runs PCR based test.
 * @param {Matrix} X - matrix with predictors.
 * @param {Matrix} Y - matrix with responses.
 * @param {ncomp} ncomp - number of components for PV-set creation.
 * @param {boolean} scale - to scale or not to scale calibration set when creating PV-set.
 * @param {JSON} cv - cross-validation resampling parameters.
 * @param {string} fileSuffix - suffix for reference files.
 */
function runPCRTest(X, Y, ncomp, scale, cv, fileSuffix) {

   // compute results for global and local CV-scopes
   const m = pcrfit(X, Y, ncomp, center, scale);

   const [Xpvg, Dg] = pcvpcr(X, Y, m, ncomp, cv, 'global');
   const [Xpvl, Dl] = pcvpcr(X, Y, m, ncomp, cv, 'local');

   const rg = pcrpredict(m, Xpvg, Y);
   const rl = pcrpredict(m, Xpvl, Y);

   // read reference values
   const tDg = readOutcome(referenceDir + 'pcvpcr/' + 'Dg' + fileSuffix);
   const tDl = readOutcome(referenceDir + 'pcvpcr/' + 'Dl' + fileSuffix);
   const tYpredg = readOutcome(referenceDir + 'pcvpcr/' + 'Ypvg' + fileSuffix);
   const tYpredl = readOutcome(referenceDir + 'pcvpcr/' + 'Ypvl' + fileSuffix);

   // tests for correct type
   expect(ismatrix(tDg) && ismatrix(Dg)).to.be.true;
   expect(ismatrix(tDl) && ismatrix(Dl)).to.be.true;
   expect(ismatrix(tYpredg) && ismatrix(rg.Ypred)).to.be.true;
   expect(ismatrix(tYpredl) && ismatrix(rl.Ypred)).to.be.true;

   // compute maximum absolute difference between computed and reference values
   const dgMaxDiff = max(tDg.subtract(Dg).apply(x => Math.abs(x), 0).v)
   const dlMaxDiff = max(tDl.subtract(Dl).apply(x => Math.abs(x), 0).v)
   const ygMaxDiff = max(tYpredg.subtract(rg.Ypred).apply(x => Math.abs(x), 0).v)
   const ylMaxDiff = max(tYpredl.subtract(rl.Ypred).apply(x => Math.abs(x), 0).v)

   // test that the difference does not exceed tolerance limits
   expect(dgMaxDiff < tolD).to.be.true;
   expect(dlMaxDiff < tolD).to.be.true;
   expect(ygMaxDiff < tolY).to.be.true;
   expect(ylMaxDiff < tolY).to.be.true;
}

/**
 * Runs PLS based test.
 * @param {Matrix} X - matrix with predictors.
 * @param {Matrix} Y - matrix with responses.
 * @param {ncomp} ncomp - number of components for PV-set creation.
 * @param {boolean} scale - to scale or not to scale calibration set when creating PV-set.
 * @param {JSON} cv - cross-validation resampling parameters.
 * @param {string} fileSuffix - suffix for reference files.
 */
function runPLSTest(X, Y, ncomp, scale, cv, fileSuffix) {

   // compute results for global and local CV-scopes
   const m = plsfit(X, Y, ncomp, center, scale);

   const [Xpvg, Dg] = pcvpls(X, Y, m, ncomp, cv, 'global');
   const [Xpvl, Dl] = pcvpls(X, Y, m, ncomp, cv, 'local');

   const rg = plspredict(m, Xpvg, Y);
   const rl = plspredict(m, Xpvl, Y);

   // read reference values
   const tDg = readOutcome(referenceDir + 'pcvpls/' + 'Dg' + fileSuffix);
   const tDl = readOutcome(referenceDir + 'pcvpls/' + 'Dl' + fileSuffix);
   const tYpredg = readOutcome(referenceDir + 'pcvpls/' + 'Ypvg' + fileSuffix);
   const tYpredl = readOutcome(referenceDir + 'pcvpls/' + 'Ypvl' + fileSuffix);

   // tests for correct type
   expect(ismatrix(tDg) && ismatrix(Dg)).to.be.true;
   expect(ismatrix(tDl) && ismatrix(Dl)).to.be.true;
   expect(ismatrix(tYpredg) && ismatrix(rg.Ypred)).to.be.true;
   expect(ismatrix(tYpredl) && ismatrix(rl.Ypred)).to.be.true;

   // compute maximum absolute difference between computed and reference values
   const dgMaxDiff = max(tDg.subtract(Dg).apply(x => Math.abs(x), 0).v)
   const dlMaxDiff = max(tDl.subtract(Dl).apply(x => Math.abs(x), 0).v)
   const ygMaxDiff = max(tYpredg.subtract(rg.Ypred).apply(x => Math.abs(x), 0).v)
   const ylMaxDiff = max(tYpredl.subtract(rl.Ypred).apply(x => Math.abs(x), 0).v)

   // test that the difference does not exceed tolerance limits
   expect(dgMaxDiff < tolD).to.be.true;
   expect(dlMaxDiff < tolD).to.be.true;
   expect(ygMaxDiff < tolY).to.be.true;
   expect(ylMaxDiff < tolY).to.be.true;
}

/**
 * Runs test for various combination of PCV algorithm parameters.
 * @param {function} fruntest - test function.
 */
function runAllTests(fruntest) {

   // cv -> string
   const getcvstr = (cv) => cv.nseg ? cv.type + cv.nseg : cv.type;

   // combination of parameters to test
   const ncomp_cases = [1, 10, 20, 30];
   const cv_cases = [{type: "ven", nseg: 4}, {type: "ven", nseg: 10}, {type: "loo"}];
   const scale_cases = [true, false];

   // const ncomp_cases = [10];
   // const cv_cases = [{type: "ven", nseg: 4}];
   // const scale_cases = [true];

   // load test data (Corn)
   const D = readOutcome(referenceDir + 'data/corn.csv');
   const [X, Y] = splitregdata(D);

   for (let ic in ncomp_cases) {
      for (let icv in cv_cases) {
         for (let isc in scale_cases) {

            const ncomp = ncomp_cases[ic];
            const scale = scale_cases[isc];
            const cv = cv_cases[icv];
            const fileSuffix = "-" + ncomp + "-" + scale.toString().toUpperCase() + "-" + getcvstr(cv) + ".csv";
            fruntest(X, Y, ncomp, scale, cv, fileSuffix);
         }
      }
   }
}


describe('Reference tests for PCV methods', function () {

   it ('tests for method "pcvpls"', function () {
      if (!checkTestDir('pcvpls')) this.skip();
      runAllTests(runPLSTest);
   }).timeout(1000000000);

   it ('tests for method "pcvpcr"', function () {
      if (!checkTestDir('pcvpcr')) this.skip();
      runAllTests(runPCRTest);
   }).timeout(1000000000);

   it ('tests for method "pcvpca"', function () {
      if (!checkTestDir('pcvpca')) this.skip();
      runAllTests(runPCATest);
   }).timeout(1000000000);

});
