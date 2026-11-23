% Build the toolbox file (.mltbx) for File Exchange.
% Run from the MATLAB folder of the repository (works in MATLAB Online too, R2023a or newer).

version = "1.2.0";
root = string(pwd);

% the same GUID as in pcv.prj, so the new file updates the installed toolbox
opts = matlab.addons.toolbox.ToolboxOptions(root, "29122b16-1d7b-4a42-aa54-b66a9654c313");
opts.ToolboxName = "Procrustes cross-validation";
opts.ToolboxVersion = version;

% package name used by File Exchange ("mpminstall procrustes_cv"), property exists from R2026b
if isprop(opts, "PackageName")
   opts.PackageName = "procrustes_cv";
end

opts.AuthorName = "Sergey Kucheryavskiy";
opts.AuthorEmail = "svkucheryavski@gmail.com";
opts.Summary = "Set of methods for generating Procrustes validation sets for PCA/SIMCA, PCR and PLS models.";
% R2026b replaced Description with Readme
if isprop(opts, "Readme")
   opts.Readme = fullfile(root, "README.md");
else
   opts.Description = "Procrustes cross-validation is a new approach for validation of chemometric models. " + ...
      "It makes possible to generate a new dataset, named ""PV-set"" and use it for validation of models " + ...
      "in the same way as with an independent validation set. The current implementation supports " + ...
      "Principal Component Analysis (method ""pcvpca()""), Principal Component Regression (method ""pcvpcr()"") " + ...
      "and Partial Least Squares (method ""pcvpls()"") models." + newline + newline + ...
      "Check Getting started documentation and see more details here: https://github.com/svkucheryavski/pcv.";
end

% files included in the toolbox (tests and project files are left out)
opts.ToolboxFiles = fullfile(root, ["README.md", "corn.mat", "crossval.m", "demo.m", "doc", ...
   "getxpvorth.m", "pcv.m", "pcv3wbu.m", "pcv3wfu.m", "pcv3whsv.m", "pcv3wtools.m", "pcvpca.m", ...
   "pcvpcr.m", "pcvpls.m", "pcvreg.m", "simpls.m"]);
opts.ToolboxGettingStartedGuide = fullfile(root, "doc", "GettingStarted.mlx");
opts.ToolboxImageFile = fullfile(root, "..", "logo.png");

% File Exchange requires the BSD license in uploaded packages (property exists from R2026b)
if isprop(opts, "UseLicenseBSD")
   opts.UseLicenseBSD = true;
else
   warning("This MATLAB release cannot add the BSD license required by File Exchange; use R2026b or newer.");
end

% the code uses "..." strings, available from R2017a
opts.MinimumMatlabRelease = "R2017a";

opts.SupportedPlatforms.Win64 = true;
opts.SupportedPlatforms.Glnxa64 = true;
opts.SupportedPlatforms.Mac = true;   % Intel and Apple silicon
opts.SupportedPlatforms.MatlabOnline = true;

opts.OutputFile = fullfile(root, "Procrustes cross-validation.mltbx");
matlab.addons.toolbox.packageToolbox(opts);
fprintf("Toolbox %s saved to %s\n", version, opts.OutputFile);
