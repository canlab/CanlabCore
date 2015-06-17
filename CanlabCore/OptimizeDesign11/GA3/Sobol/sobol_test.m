%% SOBOL_TEST runs the Sobol tests.
%
timestamp;
fprintf ( 1, '\n' );
fprintf ( 1, 'SOBOL_TEST\n' );
fprintf ( 1, '  Test the MATLAB SOBOL routines.\n' );

sobol_test01;
sobol_test02;
sobol_test03;
sobol_test04;
sobol_test05;

fprintf ( 1, '\n' );
fprintf ( 1, 'SOBOL_TEST\n' );
fprintf ( 1, '  Normal end of execution.\n' );

fprintf ( 1, '\n' );
timestamp;
