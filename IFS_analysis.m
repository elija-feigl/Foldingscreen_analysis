%% execute this to select a folder and analyze a folding screen
[pname] = uigetdir(pwd, 'Select an Directory with intial foling screen');
sigma_integrate_band = 1.0;

[name] = step1_make_profiles(pname);
step2_anotate_profiles(name, pname, sigma_integrate_band);
step3_analyse_profiles(name, pname);
