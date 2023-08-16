/*****************************************************************
 *  Speed tests for PCV methods                                  *
 *****************************************************************/




// import dependencies
import * as fs from 'fs';
import { Matrix } from 'mdatools/arrays';
import { pcafit, pcapredict, pcrfit, pcrpredict, plsfit, plspredict, splitregdata } from 'mdatools/models';

// import methods to test
import { pcvpca, pcvpcr, pcvpls } from '../src/index.js';


function measure(f, msg = '') {
   const start = Date.now();
   const out = f()
   console.log(msg + ' - done in: ' + (Date.now() - start) / 1000 + ' s.')
   return out;
}

const A = 30;
const cv = {'type': 'ven', nseg: 10}
const D = Matrix.parseCSV(fs.readFileSync('../.tests/data/corn.csv', 'utf8')).values;
const [X, Y] = splitregdata(D);

console.log(`\nSpeed tests for PCV methods with ${X.nrows} x ${X.ncols} matrix (A = ${A}, nseg = ${cv.nseg}):`);
console.log('---------------------------------------------------------------------');

console.log('');
const mpca = measure( () => pcafit(X, A), 'pcafit');
measure( () => pcapredict(mpca, X), 'pcapredict');
measure( () => pcvpca(X, mpca, A, cv), 'pcvpca (fast)');
measure( () => pcvpca(X, mpca, A, cv, 'global', true), 'pcvpca (precise)');

console.log('');
const mpcr = measure( () => pcrfit(X, Y, A), 'pcrfit');
measure( () => pcrpredict(mpcr, X), 'pcrpredict');
measure( () => pcvpcr(X, Y, mpcr, A, cv), 'pcvpcr (fast)');
measure( () => pcvpcr(X, Y, mpcr, A, cv, 'global', true), 'pcvpcr (precise)');

console.log('');
const mpls = measure( () => plsfit(X, Y, A), 'plsfit');
measure( () => plspredict(mpls, X, Y), 'plspredict');
measure( () => pcvpls(X, Y, mpls, A, cv), 'pcvpls (fast)');
measure( () => pcvpls(X, Y, mpls, A, cv, 'global', true), 'pcvpls (precise)');

console.log('');
