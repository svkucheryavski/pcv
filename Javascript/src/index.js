/************************************************************/
/*    Methods for Procrustes Cross Validation               */
/************************************************************/


// import methods to test
import { _dot, ismatrix, crossprod, tcrossprod, Vector, Index, Matrix, isindex} from 'mdatools/arrays';
import { ssq, norm2, max } from 'mdatools/stat';
import { rsvd } from 'mdatools/decomp';
import { scale as prep_scale, unscale as prep_unscale } from 'mdatools/prep';


/**
 * Returns settings for RSVD algorithm depending on criterion (fast or precise)
 * @param {boolean} precise - shall precise version be selected or not.
 * @returns {Array} array with settings for its, pa, and pb
 */
function get_rsvd_settings(precise) {
   if (precise) return [1.2, 15, 4] // [pa, pb, its]
   return [1, 15, 3];
}


/**
 * Compute PV-set using PLS.
 *
 * @param {Matrix} data - matrix with calibration set (first column with y the rest is X).
 * @param {JSON} m - object with global PCR model.
 * @param {number} ncomp - number of components to use.
 * @param {JSON} cv - object with type of splits ('ven', 'rand', 'loo') and number of segments.
 * @param {string} [cvscope="global"] - defines if autoscale must be done globally ('global') or locally ('local').
 * @param {boolean} [precise=false] - shall precise version be selected or not.
 *
 * @return {Matrix} matrix with PV-set.
 *
 */
export function pcvpls(X, Y, m, ncomp, cv, cvscope, precise) {

   // assemble JSON with functions
   const funlist = {

      // adjust a global PCR model
      getglobalmodel: function(m, cvncomp) {
         const compind = Index.seq(1, cvncomp);
         const R = m.ncomp === cvncomp ? m.R : m.R.subset([], compind);
         const P = m.ncomp === cvncomp ? m.P : m.P.subset([], compind);
         const C = m.ncomp === cvncomp ? m.C : m.C.subset([], compind);

         return {ncomp: cvncomp, R: R, P: P, C: C, PRM: Matrix.eye(X.ncols, X.ncols).subtract(tcrossprod(m.R, m.P))}
      },

      // creates a local PCR model
      getlocalmodel: function(Xc, Yc, m) {

         // get loadings for local model
         const mk = simpls(Xc, Yc, m.ncomp, precise);

         // correct direction of loadings for local model and get scores
         const a = getdirections(m.R, mk.R);
         const A = Matrix.diagm(a, a.length);
         const Pka = mk.P.dot(A);
         const Rka = mk.R.dot(A);
         const Cka = mk.C.dot(A);

         return {P: Pka, C: Cka, R: Rka}
      },

      // compute the explained part of Xpv
      getxpv: function(m, mk, Xk) {
         const Tk = Xk.dot(mk.R);
         const dk = mk.C.divide(m.C);

         const Tpvk = Matrix.zeros(Tk.nrows, Tk.ncols);
         for (let a = 1; a <= m.ncomp; a++) {
            const tpva = Tpvk.getcolref(a);
            const ta = Tk.getcolref(a);
            for (let i = 0; i < Tk.nrows; i++) {
               tpva[i] = ta[i] * dk.v[a - 1];
            }
         }

         const Xpvk = tcrossprod(Tpvk, m.P);
         return [Xpvk, dk];
      },

      // compute orthogonal distances to local model
      getqk: function(Xk, mk) {
         const PRM = Matrix.eye(Xk.ncols, Xk.ncols).subtract(tcrossprod(mk.R, mk.P));
         const Ek = Xk.dot(PRM);
         return Ek.apply(ssq, 1);
      }
   }

   return pcvreg(X, Y, m, ncomp, cv, funlist, cvscope);
}


/**
 * Compute PV-set using Principal Component Regression.
 *
 * @param {Matrix} X - matrix with calibration set predictors.
 * @param {Matrix} Y - matrix with calibration set responses (1 column).
 * @param {JSON} m - object with global PCR model.
 * @param {JSON} cv - object with type of splits ('ven', 'rand', 'loo') and number of segments.
 * @param {string} [cvscope="global"] - defines if autoscale must be done globally ('global') or locally ('local').
 * @param {boolean} [precise=false] - shall precise version be selected or not.
 *
 * @return {Array} two matrices: Xpv with PV-set and D with scalars.
 *
 * @return {Array} array with two matrices: D with scalars and Xpv with PV-set.
 *
 */
export function pcvpcr(X, Y, m, ncomp, cv, cvscope, precise) {

   const [pa, pb, its] = get_rsvd_settings(precise);

   // assemble JSON with functions
   const funlist = {

      // adjust a global PCR model
      getglobalmodel: function(m, cvncomp) {
         const P = m.ncomp === cvncomp ? m.P : m.P.subset([], Index.seq(1, cvncomp));
         const C = m.ncomp === cvncomp ? m.C : m.C.subset([], Index.seq(1, cvncomp));

         return {ncomp: cvncomp, P: P, C: C, PRM: Matrix.eye(X.ncols, X.ncols).subtract(tcrossprod(P, P))}
      },

      // creates a local PCR model
      getlocalmodel: function(Xc, Yc, m) {

         // get loadings for local model
         const Pk = rsvd(Xc, m.ncomp, pa, pb, its).V;

         // correct direction of loadings for local model and get scores
         const a = getdirections(m.P, Pk);
         const Pka = Pk.dot(Matrix.diagm(a, a.length));
         const Tc = Xc.dot(Pka);

         // compute Y-loadings for local model
         const Ck = Matrix.zeros(1, m.ncomp);
         for (let a = 1; a <= m.ncomp; a++) {
            const ta = Tc.getcolref(a);
            const eigena = ssq(ta);
            Ck.v[a - 1] = _dot(Yc.v, ta, 1, Xc.nrows, Xc.nrows, 1) / eigena;
         }

         return {P: Pka, C: Ck}
      },

      // compute the explained part of Xpv
      getxpv: function(m, mk, Xk) {
         const Tk = Xk.dot(mk.P);
         const dk = mk.C.divide(m.C);

         const Tpvk = Matrix.zeros(Tk.nrows, Tk.ncols);
         for (let a = 1; a <= m.ncomp; a++) {
            const tpva = Tpvk.getcolref(a);
            const ta = Tk.getcolref(a);
            for (let i = 0; i < Tk.nrows; i++) {
               tpva[i] = ta[i] * dk.v[a - 1];
            }
         }

         const Xpvk = tcrossprod(Tpvk, m.P);
         return [Xpvk, dk];
      },

      // compute orthogonal distances to local model
      getqk: function(Xk, mk) {
         const PI = Matrix.eye(Xk.ncols, Xk.ncols).subtract(tcrossprod(mk.P, mk.P));
         const Ek = Xk.dot(PI);
         return Ek.apply(ssq, 1);
      }
   }

   return pcvreg(X, Y, m, ncomp, cv, funlist, cvscope);
}


/**
 * Compute PV-set for regression models.
 *
 * @param {Matrix} X - matrix with predictors from calibration set.
 * @param {Matrix} Y - matrix (1 column) with responses from calibration set.
 * @param {JSON} m - object with global model.
 * @param {JSON} cv - object with type of splits ('ven', 'rand', 'loo') and number of segments.
 * @param {JSON} funlist - list of functions to compute different parts of Xpv.
 * @param {string} [cvscope="global"] - defines if autoscale must be done globally ('global') or locally ('local').
 *
 * @return {Array} two matrices: Xpv with PV-set and D with scalars.
 *
 */
function pcvreg(X, Y, m, ncomp, cv, funlist, cvscope) {

   if (!cvscope) {
      cvscope = 'global';
   }

   if (!ismatrix(X) || X.nrows < 4 || X.ncols < 2 ) {
      throw new Error('pcvreg: parameter "X" must be a matrix with at least 4 rows and 2 columns.')
   }

   if (!ismatrix(Y) || Y.nrows !== X.nrows ) {
      throw new Error('pcvreg: parameter "Y" must be a matrix with the same number of rows as "X".')
   }

   if (cv === null || typeof(cv) !== 'object') {
      throw new Error('pcvreg: parameter "cv" must be specified correctly.')
   }

   if (ncomp === undefined || ncomp === null) {
      ncomp = m.ncomp;
   }

   if (ncomp < 1 || ncomp > X.nrows) {
      throw new Error('pcvreg: wrong value for parameter "ncomp".')
   }

   if (ncomp > m.ncomp) {
      throw new Error(`pcvreg: value for parameter "ncomp" can not exceed number of components in global model "m" (${m.ncomp}).`)
   }

   const nobs = X.nrows;
   const nvar = X.ncols;

   // get cross-validation parameters and adjusted number of components
   const [cvncomp, cvind, cvnseg] = getcvparams(nobs, nvar, cv, ncomp, Y.getcolumn(1));

   // create object with adjusted global model
   const mg = funlist.getglobalmodel(m, cvncomp);

   // autoscale the calibration set globally
   const Xp = prep_scale(X, m.mX, m.sX);
   const Yp = prep_scale(Y, m.mY, m.sY);

   // prepare empty matrices for pseudo-validation set and scalars
   const Xpv = Matrix.zeros(nobs, nvar);
   const D = Matrix.zeros(cvncomp, cvnseg);

   // loop for computing the PV set
   for (let k = 1; k <= cvnseg; k++) {

      // get row indices for local calibration and validation sets
      const [indc, indk] = cv2obs(cvind, k);

      // get the sets
      let Xc = Xp.subset(indc, []);
      let Xk = Xp.subset(indk, []);
      let Yc = Yp.subset(indc, []);

      // center and scale the sets locally if needed
      if (cvscope == 'local') {
         let mX, sX;
         [Xc, mX, sX] = prep_scale(Xc, m.center, m.scale, true);
         Yc = prep_scale(Yc, m.center, m.scale);
         Xk = prep_scale(Xk, mX, sX);
      }

      const mk = funlist.getlocalmodel(Xc, Yc, mg);
      const [Xpvhat, dk] = funlist.getxpv(mg, mk, Xk);
      const Xpvorth = getxpvorth(funlist.getqk(Xk, mk), Xk, mg.PRM);

      Xpv.replace(Xpvhat.add(Xpvorth), indk, []);
      D.v.set(dk.v, (k - 1) * cvncomp)
   }

   // uncenter and unscale the values and return
   return [prep_unscale(Xpv, m.mX, m.sX), D.t()];
}


/**
 * Compute PV-set using PCA.
 *
 * @param {Matrix} X - matrix with calibration set.
 * @param {JSON} m - object with global PCA model.
 * @param {number} ncomp - number of components to use.
 * @param {JSON} cv - object with type of splits ('ven', 'rand', 'loo') and number of segments.
 * @param {string} [cvscope="global"] - defines if autoscale must be done globally ('global') or locally ('local').
 * @param {boolean} [precise=false] - shall precise version be selected or not.
 *
 * @return {Matrix} matrix with PV-set.
 *
 */
export function pcvpca(X, m, ncomp, cv, cvscope, precise) {

   if (!cvscope) {
      cvscope = 'global';
   }

   if (!ismatrix(X) || X.nrows < 4 || X.ncols < 2 ) {
      throw new Error('pcvpca: parameter "X" must be a matrix with at least 4 rows and 2 columns.')
   }

   if (cv === null || typeof(cv) !== 'object') {
      throw new Error('pcvpca: parameter "cv" must be specified correctly.')
   }

   if (ncomp === undefined || ncomp === null) {
      ncomp = m.ncomp;
   }

   if (ncomp < 1 || ncomp > X.nrows) {
      throw new Error('pcvpca: wrong value for parameter "ncomp".')
   }

   if (ncomp > m.ncomp) {
      throw new Error(`pcvpca: value for parameter "ncomp" can not exceed number of components in global model "m" (${m.ncomp}).`)
   }

   const [pa, pb, its] = get_rsvd_settings(precise);
   const nobs = X.nrows;
   const nvar = X.ncols;

   // get cross-validation parameters and adjusted number of components
   const [cvncomp, cvind, cvnseg] = getcvparams(nobs, nvar, cv, ncomp, m);

   // autoscale the calibration set
   const Xp = prep_scale(X, m.mX, m.sX);

   // create a global model
   const P = m.P.subset([], Index.seq(1, cvncomp));
   const PRM = Matrix.eye(nvar).subtract(tcrossprod(P));

   // prepare empty matrix for pseudo-validation set
   const Xpv = Matrix.zeros(nobs, nvar);

   // loop for computing the PV set
   for (let k = 1; k <= cvnseg; k++) {

      // get row indices for local calibration and validation sets
      const [indc, indk] = cv2obs(cvind, k);

      // get the sets
      let Xc = Xp.subset(indc, []);
      let Xk = Xp.subset(indk, []);

      // center and scale the sets locally if needed
      if (cvscope == 'local') {
         let mX, sX;
         [Xc, mX, sX] = prep_scale(Xc, m.center, m.scale, true);
         Xk = prep_scale(Xk, mX, sX);
      }

      // get loadings for local model and rotation matrix between global and local models
      const Pk = rsvd(Xc, cvncomp, pa, pb, its).V;

      // correct direction of loadings for local model
      const a = getdirections(P, Pk);
      const Pka = Pk.dot(Matrix.diagm(a, nvar));

      // compute explained part of Xpvk
      const Tk = Xk.dot(Pka);
      const Xpvk_hat = tcrossprod(Tk, P);


      // if decomposition is full return the explained part only
      if (P.nrows > P.ncols) {
         const Ek = Xk.subtract(tcrossprod(Tk, Pka));
         const qkn = Ek.apply(norm2, 1);
         const Xpvk_orth = getxpvorth(qkn, Xk, PRM);
         Xpv.replace(Xpvk_hat.add(Xpvk_orth), indk, []);
      } else {
         Xpv.replace(Xpvk_hat, indk, []);
      }
   }

   // ucenter and unscale the values and return
   return prep_unscale(Xpv, m.mX, m.sX);
}


/**
 * Compute orthogonal part of PV-set.
 *
 * @param {Vector} qkn - vector with orthogonal distances for k-th segment (square root of qk).
 * @param {Matrix} Xk - matrix with local validation set for k-th segment.
 * @param {Matrix} PRM - projection matrix for global model.
 *
 * @returns {Matrix} matrix with orthogonal part values.
 *
 */
export function getxpvorth(qkn, Xk, PRM) {


   const nobj = qkn.length;

   // project Xk to random vectors
   const Z = Matrix.rand(nobj, nobj, -1, 1)
   const X1 = Z.dot(Xk);

   // normalize columns of X1 to norm2
   for (let c = 1; c <= X1.ncols; c++) {
      const x = X1.getcolref(c);
      const nx = norm2(x);
      for (let r = 1; r <= nobj; r++) {
         x[r - 1] = nx == 0 ? x[r - 1] : x[r - 1] / nx;
      }
   }

   // compute X2 and scale its rows to unit length
   const X2 = X1.dot(PRM);

   // normalize rows of X2 to sqrt(qk/norm2)
   const nr = X2.apply(norm2, 1);
   for (let c = 1; c <= X1.ncols; c++) {
      const x = X2.getcolref(c);
      for (let r = 1; r <= nobj; r++) {
         x[r - 1] = nr.v[r - 1] == 0 ? x[r - 1] : x[r - 1] * qkn.v[r - 1] / nr.v[r - 1];
      }
   }

   return X2;
}


/**
 * Return two vector with observation indices (for local training and for local validation subsets).
 *
 * @param {Index} cv - vector with CV segments indices (e.g. from 'crossval()' method).
 * @param {number} k - which segment to use for validation.
 *
 * @return {Array} - array with two Int32Array with indices (first for "cal" and second for "val").
 *
 */
export function cv2obs(cv, k) {

   if (!isindex(cv)) {
      throw Error('cv2obs: parameter "cv" must be instance of Index.');
   }

   return [
      cv.which(v => v !== k),
      cv.which(v => v === k)
   ];
}


/**
 * Generate indices for cross-validation.
 *
 * @param {string} type - cross-validation type ("loo", "ven", "rand").
 * @param {number} nseg - number of segments.
 * @param {number} nobs - number of observations.
 *
 * @return {Index} vector of length 'nobs' with segment indices for each observation.
 *
 */
export function crossval(type, nseg, nobs, yref) {

   if (type == "loo") {
      return Index.seq(1, nobs);
   }

   if (!yref || yref.length !== nobs) {
      yref = Vector.seq(1, nobs);
   }

   if (type == "ven") {
      const nrep = Math.ceil(nobs/nseg);
      const cvInd = new Index(Index.seq(1, nseg).rep(nrep).v.subarray(0, nobs));
      const yInd = yref.sortind();
      return cvInd.subset(yInd.sortind());
   }

   if (type == "rand") {
      return crossval("ven", nseg, nobs).shuffle();
   }
}


/**
 * Check if direction of vectors from 'X' and 'Xk' are same or opposite.
 *
 * @param {Matrix} X - matrix with vectors for global model (e.g. loadings).
 * @param {Matrix} Xk - matrix with vectors for local model (e.g. loadings).
 *
 * @returns {Vector} vector with '1' if directions are the same and '-1' if opposite.
 *
 */
export function getdirections(X, Xk) {

   const a = Vector.zeros(X.ncols);

   for (let c = 1; c <= X.ncols; c++) {
      const xc = X.getcolref(c);
      const xck = Xk.getcolref(c);
      const cnorm = norm2(xc);
      const cnormk = norm2(xck);

      let s = 0;
      for (let r = 0; r < xc.length; r++) {
         s += xc[r] / cnorm * xck[r] / cnormk;
      }

      a.v[c - 1] = (Math.acos(s) < (Math.PI / 2.0)) * 2 - 1;
   }

   return a;
}


/**
 * Computes cross-validation parameters for a model.
 *
 * @param {number} nobs - number of observations.
 * @param {number} nvar - number of predictors.
 * @param {number} ncomp - number of components specified by user.
 * @param {Object} cv - JSON wiht user defined cross-validation parameters.
 *
 * @returns {Array} array with four elements: adjusted number of components, type, vector with
 * indices and number of segments in cross-validation.
 *
 */
export function getcvparams(nobs, nvar, cv, ncomp, resp) {

   const nseg = cv === undefined ? 4 : cv.nseg;
   const cvType = cv === undefined ? "ven" : cv.type;
   const cvInd = crossval(cvType, nseg, nobs, resp);
   const cvNSeg = max(cvInd);

   // correct maximum number of components
   const maxNcomp = Math.min(ncomp, nobs - Math.ceil(nobs / nseg) - 1 , nvar);
   if (ncomp > maxNcomp) {
      ncomp = maxNcomp;
   }

   return [ncomp, cvInd, cvNSeg]
}


/**
 * Implementation of SIMPLS algorithm.
 *
 * @param {Matrix} X - matrix with predictors.
 * @param {Matrix} Y - matrix with responses.
 * @param {number} ncomp - number of components.
 *
 * @return {JSON} JSON with decomposition results.
 *
 */
export function simpls(X, Y, ncomp, precise) {

   const [pa, pb, its] = get_rsvd_settings(precise);

   const nobj  = X.nrows;
   const npred = X.ncols;
   const nresp = Y.ncols;

   // initial estimation
   let S = crossprod(X, Y);
   let M = crossprod(X);

   // prepare space for results
   const C = Matrix.zeros(nresp, ncomp)
   const R = Matrix.zeros(npred, ncomp)
   const V = Matrix.zeros(npred, ncomp)
   const P = Matrix.zeros(npred, ncomp)
   const T = Matrix.zeros(nobj, ncomp)
   const U = Matrix.zeros(nobj, ncomp)

   const xeigenvals = Vector.zeros(ncomp);
   const yeigenvals = Vector.zeros(ncomp);

   // loop for each components
   for (let a = 1; a <= ncomp; a++) {

      let r = rsvd(S, 1, pa, pb, its).U
      let t = X.dot(r);

      const tnorm = norm2(t.v);
      t = t.apply(v => v / tnorm, 0);
      r = r.apply(v => v / tnorm, 0);

      const p = crossprod(X, t);
      const c = crossprod(Y, t);
      let u = Y.dot(c);
      let v = p.copy()

      if (a > 1) {
         v = v.subtract(V.dot(crossprod(V, p)));
         u = u.subtract(T.dot(crossprod(T, u)));
      }

      const vnorm = norm2(v.v);
      v = v.apply(x => x / vnorm, 0);

      R.v.set(r.v, (a - 1) * npred);
      V.v.set(v.v, (a - 1) * npred);
      P.v.set(p.v, (a - 1) * npred);
      T.v.set(t.v, (a - 1) * nobj);
      U.v.set(u.v, (a - 1) * nobj);
      C.v.set(c.v, (a - 1) * nresp);


      xeigenvals.v[a - 1] = ssq(t.v) / (nobj - 1);
      yeigenvals.v[a - 1] = ssq(u.v) / (nobj - 1);

      M = M.subtract(tcrossprod(p))
      S = S.subtract(v.dot(crossprod(v, S)));
   }

   return {
      R: R, P: P, T: T, C: C, U: U,
      xeigenvals: xeigenvals,
      yeigenvals: yeigenvals
   };
}
