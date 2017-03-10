function isitinstalled=wekaInstalled
s=javaclasspath;

isitinstalled=~isempty(strfind(s,'weka.jar'));
