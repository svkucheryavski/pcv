% Compute vector with segment indices for cross-validation
%
% Arguments
% ---------
% CV        number of segments or a cell array (see description)
% nObj      number of objects (rows) in training set
% Y         matrix or vector with response values (needed for "ven" splits)
%
% The method computes vector with segment number for each row. The 'CV' 
% parameter can be one of the follows:
% 
% 'CV = 1' - full cross-validation (leave one out)
% 'CV = n' - random cross-validation with n segments
% 'CV = {"rand", n}' - random cross-validation with n segments
% 'CV = {"ven", n}' - systematic cross-validation with n segments
% 'CV = {"loo"}' - full cross-validation (leave one out)
%
function ind = crossval(CV, nObj, Y) 
   
   % if user already provided vector with values - return it
   if isa(CV, 'double') && numel(CV) == nObj
      ind = CV;
      return;
   end   

   % check and process CV parameters
   p = getcvparams(CV, nObj);

   % check number of segments
   if p{2} < 2 || p{2} > nObj
      error('Wrong value for number of segments (should be between 2 and number of objects).')
   end

   % leave-one-out
   if p{1} == "loo"
      ind = (1:nObj)';
      return;
   end

   nSeg = p{2};
   nMax = ceil(nObj / nSeg); 
   ind = repmat(1:nSeg, 1, nMax);
   ind = ind(1:nObj)';
   
   % systematic
   if p{1} == "ven"
      if nargin < 3
         return;
      end   
      [~, rowInd] = sort(Y(:, 1));
      ind = ind(rowInd);
      return;
   end

   % random
   if p{1} == "rand"
      ind = ind(randperm(numel(ind)));
      return;
   end   

   error("Something went wrong.");
end


function params = getcvparams(CV, nObj)

   % a number - random split
   if isa(CV, 'double')
      params = {"rand", (CV == 1) * nObj + (CV > 1) * CV};
      return;
   end

   type = CV{1};

   % leave one out
   if type == "loo"
      params = {"loo", nObj};
      return;
   end   

   % venetian blinds
   if type == "ven" || type == "rand"
      params = CV;
      return;
   end

   error("Wrong value for 'CV' parameter.");
end