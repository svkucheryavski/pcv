/*****************************************************************
 *  Speed tests for PCV methods                                  *
 *****************************************************************/

// import dependencies
import { Matrix } from 'mdatools/arrays';
import { pcafit, pcapredict, pcrfit, pcrpredict, plsfit, plspredict } from 'mdatools/models';

// import methods to test
import { pcvpca, pcvpcr, pcvpls } from '../src/index.js';


function measure(f, msg = '') {
   const start = Date.now();
   const out = f()
   console.log(msg + ' - done in: ' + (Date.now() - start) / 1000 + ' s.')
   return out;
}

const X = Matrix.rand(100, 1000);
const Y = Matrix.rand(100, 1);

console.log('\nSpeed tests for PCV methods with 100 x 1000 matrix:');
console.log('---------------------------------------------------');

const mpca = measure( () => pcafit(X, 20), 'pcafit');
measure( () => pcapredict(mpca, X), 'pcapredict');
measure( () => pcvpca(X, mpca, 20, {type: 'ven', nseg: 10}), 'pcvpca');

const mpcr = measure( () => pcrfit(X, Y, 20), 'pcrfit');
measure( () => pcrpredict(mpcr, X), 'pcrpredict');
measure( () => pcvpcr(X, Y, mpcr, 20, {type: 'ven', nseg: 10}), 'pcvpcr');

const mpls = measure( () => plsfit(X, Y, 20), 'plsfit');
measure( () => plspredict(mpls, X, Y), 'plspredict');
measure( () => pcvpls(X, Y, mpls, 20, {type: 'ven', nseg: 10}), 'pcvpls');

console.log('');
