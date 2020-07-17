%% RCP data 
fid = fopen('ch4.txt');
ch4data = fscanf(fid, '%f', [736 4]);

CH425data = ch4data(:,1);
CH445data = ch4data(:,2);
CH460data = ch4data(:,3);
CH485data = ch4data(:,4);

CO285data = [2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.89E+02	2.89E+02	2.89E+02	2.90E+02	2.90E+02	2.91E+02	2.91E+02	2.92E+02	2.92E+02	2.93E+02	2.93E+02	2.93E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.96E+02	2.96E+02	2.96E+02	2.96E+02	2.97E+02	2.97E+02	2.98E+02	2.98E+02	2.99E+02	2.99E+02	2.99E+02	3.00E+02	3.00E+02	3.00E+02	3.01E+02	3.01E+02	3.01E+02	3.02E+02	3.02E+02	3.02E+02	3.03E+02	3.03E+02	3.03E+02	3.04E+02	3.04E+02	3.05E+02	3.05E+02	3.05E+02	3.06E+02	3.06E+02	3.07E+02	3.07E+02	3.08E+02	3.08E+02	3.09E+02	3.09E+02	3.09E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.11E+02	3.11E+02	3.11E+02	3.12E+02	3.12E+02	3.12E+02	3.13E+02	3.14E+02	3.14E+02	3.15E+02	3.16E+02	3.16E+02	3.17E+02	3.18E+02	3.18E+02	3.19E+02	3.20E+02	3.21E+02	3.22E+02	3.23E+02	3.24E+02	3.25E+02	3.26E+02	3.27E+02	3.29E+02	3.30E+02	3.31E+02	3.32E+02	3.33E+02	3.35E+02	3.37E+02	3.38E+02	3.40E+02	3.41E+02	3.42E+02	3.44E+02	3.45E+02	3.47E+02	3.49E+02	3.51E+02	3.52E+02	3.54E+02	3.55E+02	3.56E+02	3.57E+02	3.58E+02	3.60E+02	3.61E+02	3.63E+02	3.65E+02	3.67E+02	3.69E+02	3.70E+02	3.73E+02	3.75E+02	3.77E+02	3.79E+02	3.81E+02	3.83E+02	3.85E+02	3.87E+02	3.89E+02	3.92E+02	3.94E+02	3.96E+02	3.99E+02	4.02E+02	4.04E+02	4.07E+02	4.10E+02	4.13E+02	4.16E+02	4.19E+02	4.22E+02	4.25E+02	4.28E+02	4.31E+02	4.35E+02	4.38E+02	4.42E+02	4.45E+02	4.49E+02	4.52E+02	4.56E+02	4.60E+02	4.64E+02	4.68E+02	4.72E+02	4.76E+02	4.81E+02	4.85E+02	4.89E+02	4.94E+02	4.99E+02	5.04E+02	5.08E+02	5.13E+02	5.19E+02	5.24E+02	5.29E+02	5.35E+02	5.41E+02	5.46E+02	5.52E+02	5.58E+02	5.64E+02	5.71E+02	5.77E+02	5.83E+02	5.90E+02	5.97E+02	6.04E+02	6.11E+02	6.18E+02	6.25E+02	6.32E+02	6.39E+02	6.47E+02	6.54E+02	6.62E+02	6.69E+02	6.77E+02	6.85E+02	6.93E+02	7.01E+02	7.09E+02	7.17E+02	7.25E+02	7.33E+02	7.42E+02	7.50E+02	7.58E+02	7.67E+02	7.75E+02	7.84E+02	7.92E+02	8.01E+02	8.10E+02	8.18E+02	8.27E+02	8.36E+02	8.45E+02	8.54E+02	8.63E+02	8.72E+02	8.81E+02	8.90E+02	8.99E+02	9.08E+02	9.17E+02	9.27E+02	9.36E+02	9.45E+02	9.54E+02	9.64E+02	9.73E+02	9.83E+02	9.92E+02	1.00E+03	1.01E+03	1.02E+03	1.03E+03	1.04E+03	1.05E+03	1.06E+03	1.07E+03	1.08E+03	1.09E+03	1.10E+03	1.11E+03	1.12E+03	1.13E+03	1.14E+03	1.15E+03	1.16E+03	1.17E+03	1.18E+03	1.19E+03	1.20E+03	1.21E+03	1.22E+03	1.23E+03	1.24E+03	1.25E+03	1.26E+03	1.27E+03	1.28E+03	1.29E+03	1.30E+03	1.31E+03	1.32E+03	1.33E+03	1.34E+03	1.35E+03	1.36E+03	1.37E+03	1.38E+03	1.39E+03	1.40E+03	1.41E+03	1.42E+03	1.43E+03	1.44E+03	1.45E+03	1.46E+03	1.47E+03	1.48E+03	1.49E+03	1.50E+03	1.51E+03	1.52E+03	1.53E+03	1.54E+03	1.55E+03	1.56E+03	1.57E+03	1.58E+03	1.58E+03	1.59E+03	1.60E+03	1.61E+03	1.62E+03	1.63E+03	1.64E+03	1.64E+03	1.65E+03	1.66E+03	1.67E+03	1.68E+03	1.68E+03	1.69E+03	1.70E+03	1.71E+03	1.71E+03	1.72E+03	1.73E+03	1.74E+03	1.74E+03	1.75E+03	1.76E+03	1.76E+03	1.77E+03	1.78E+03	1.78E+03	1.79E+03	1.79E+03	1.80E+03	1.81E+03	1.81E+03	1.82E+03	1.82E+03	1.83E+03	1.83E+03	1.84E+03	1.84E+03	1.85E+03	1.85E+03	1.86E+03	1.86E+03	1.87E+03	1.87E+03	1.88E+03	1.88E+03	1.89E+03	1.89E+03	1.89E+03	1.90E+03	1.90E+03	1.91E+03	1.91E+03	1.91E+03	1.92E+03	1.92E+03	1.92E+03	1.92E+03	1.93E+03	1.93E+03	1.93E+03	1.94E+03	1.94E+03	1.94E+03	1.94E+03	1.94E+03	1.95E+03	1.95E+03	1.95E+03	1.95E+03	1.95E+03	1.95E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03	1.96E+03];
CO260data = [2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.89E+02	2.89E+02	2.89E+02	2.90E+02	2.90E+02	2.91E+02	2.91E+02	2.92E+02	2.92E+02	2.93E+02	2.93E+02	2.93E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.96E+02	2.96E+02	2.96E+02	2.96E+02	2.97E+02	2.97E+02	2.98E+02	2.98E+02	2.99E+02	2.99E+02	2.99E+02	3.00E+02	3.00E+02	3.00E+02	3.01E+02	3.01E+02	3.01E+02	3.02E+02	3.02E+02	3.02E+02	3.03E+02	3.03E+02	3.03E+02	3.04E+02	3.04E+02	3.05E+02	3.05E+02	3.05E+02	3.06E+02	3.06E+02	3.07E+02	3.07E+02	3.08E+02	3.08E+02	3.09E+02	3.09E+02	3.09E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.11E+02	3.11E+02	3.11E+02	3.12E+02	3.12E+02	3.12E+02	3.13E+02	3.14E+02	3.14E+02	3.15E+02	3.16E+02	3.16E+02	3.17E+02	3.18E+02	3.18E+02	3.19E+02	3.20E+02	3.21E+02	3.22E+02	3.23E+02	3.24E+02	3.25E+02	3.26E+02	3.27E+02	3.29E+02	3.30E+02	3.31E+02	3.32E+02	3.33E+02	3.35E+02	3.37E+02	3.38E+02	3.40E+02	3.41E+02	3.42E+02	3.44E+02	3.45E+02	3.47E+02	3.49E+02	3.51E+02	3.52E+02	3.54E+02	3.55E+02	3.56E+02	3.57E+02	3.58E+02	3.60E+02	3.61E+02	3.63E+02	3.65E+02	3.67E+02	3.69E+02	3.70E+02	3.73E+02	3.75E+02	3.77E+02	3.79E+02	3.81E+02	3.83E+02	3.85E+02	3.87E+02	3.89E+02	3.91E+02	3.93E+02	3.95E+02	3.97E+02	3.99E+02	4.01E+02	4.03E+02	4.05E+02	4.07E+02	4.09E+02	4.11E+02	4.13E+02	4.15E+02	4.17E+02	4.19E+02	4.21E+02	4.23E+02	4.25E+02	4.27E+02	4.29E+02	4.31E+02	4.33E+02	4.35E+02	4.37E+02	4.39E+02	4.41E+02	4.44E+02	4.46E+02	4.48E+02	4.51E+02	4.53E+02	4.56E+02	4.58E+02	4.61E+02	4.63E+02	4.66E+02	4.69E+02	4.72E+02	4.75E+02	4.78E+02	4.81E+02	4.84E+02	4.87E+02	4.90E+02	4.93E+02	4.97E+02	5.00E+02	5.03E+02	5.07E+02	5.11E+02	5.14E+02	5.18E+02	5.22E+02	5.26E+02	5.29E+02	5.33E+02	5.37E+02	5.41E+02	5.46E+02	5.50E+02	5.54E+02	5.58E+02	5.63E+02	5.67E+02	5.72E+02	5.76E+02	5.81E+02	5.85E+02	5.90E+02	5.94E+02	5.99E+02	6.04E+02	6.08E+02	6.12E+02	6.17E+02	6.21E+02	6.25E+02	6.28E+02	6.32E+02	6.36E+02	6.39E+02	6.43E+02	6.46E+02	6.50E+02	6.53E+02	6.56E+02	6.60E+02	6.63E+02	6.66E+02	6.70E+02	6.73E+02	6.76E+02	6.80E+02	6.83E+02	6.86E+02	6.89E+02	6.92E+02	6.94E+02	6.97E+02	7.00E+02	7.02E+02	7.05E+02	7.07E+02	7.10E+02	7.12E+02	7.14E+02	7.16E+02	7.19E+02	7.21E+02	7.23E+02	7.24E+02	7.26E+02	7.28E+02	7.30E+02	7.31E+02	7.33E+02	7.34E+02	7.36E+02	7.37E+02	7.39E+02	7.40E+02	7.41E+02	7.42E+02	7.43E+02	7.44E+02	7.45E+02	7.46E+02	7.47E+02	7.48E+02	7.48E+02	7.49E+02	7.50E+02	7.50E+02	7.51E+02	7.51E+02	7.51E+02	7.51E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02	7.52E+02];
CO245data = [2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.89E+02	2.89E+02	2.89E+02	2.90E+02	2.90E+02	2.91E+02	2.91E+02	2.92E+02	2.92E+02	2.93E+02	2.93E+02	2.93E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.96E+02	2.96E+02	2.96E+02	2.96E+02	2.97E+02	2.97E+02	2.98E+02	2.98E+02	2.99E+02	2.99E+02	2.99E+02	3.00E+02	3.00E+02	3.00E+02	3.01E+02	3.01E+02	3.01E+02	3.02E+02	3.02E+02	3.02E+02	3.03E+02	3.03E+02	3.03E+02	3.04E+02	3.04E+02	3.05E+02	3.05E+02	3.05E+02	3.06E+02	3.06E+02	3.07E+02	3.07E+02	3.08E+02	3.08E+02	3.09E+02	3.09E+02	3.09E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.11E+02	3.11E+02	3.11E+02	3.12E+02	3.12E+02	3.12E+02	3.13E+02	3.14E+02	3.14E+02	3.15E+02	3.16E+02	3.16E+02	3.17E+02	3.18E+02	3.18E+02	3.19E+02	3.20E+02	3.21E+02	3.22E+02	3.23E+02	3.24E+02	3.25E+02	3.26E+02	3.27E+02	3.29E+02	3.30E+02	3.31E+02	3.32E+02	3.33E+02	3.35E+02	3.37E+02	3.38E+02	3.40E+02	3.41E+02	3.42E+02	3.44E+02	3.45E+02	3.47E+02	3.49E+02	3.51E+02	3.52E+02	3.54E+02	3.55E+02	3.56E+02	3.57E+02	3.58E+02	3.60E+02	3.61E+02	3.63E+02	3.65E+02	3.67E+02	3.69E+02	3.70E+02	3.73E+02	3.75E+02	3.77E+02	3.79E+02	3.81E+02	3.83E+02	3.85E+02	3.87E+02	3.89E+02	3.91E+02	3.93E+02	3.96E+02	3.98E+02	4.00E+02	4.02E+02	4.04E+02	4.07E+02	4.09E+02	4.11E+02	4.13E+02	4.16E+02	4.18E+02	4.20E+02	4.23E+02	4.25E+02	4.28E+02	4.30E+02	4.33E+02	4.35E+02	4.38E+02	4.40E+02	4.43E+02	4.45E+02	4.48E+02	4.50E+02	4.53E+02	4.56E+02	4.58E+02	4.61E+02	4.63E+02	4.66E+02	4.69E+02	4.71E+02	4.74E+02	4.76E+02	4.79E+02	4.81E+02	4.84E+02	4.87E+02	4.89E+02	4.92E+02	4.94E+02	4.96E+02	4.98E+02	5.01E+02	5.03E+02	5.05E+02	5.07E+02	5.09E+02	5.11E+02	5.13E+02	5.14E+02	5.16E+02	5.18E+02	5.19E+02	5.20E+02	5.22E+02	5.23E+02	5.24E+02	5.25E+02	5.27E+02	5.27E+02	5.28E+02	5.29E+02	5.30E+02	5.30E+02	5.31E+02	5.31E+02	5.31E+02	5.31E+02	5.31E+02	5.32E+02	5.32E+02	5.32E+02	5.32E+02	5.33E+02	5.33E+02	5.33E+02	5.34E+02	5.34E+02	5.35E+02	5.35E+02	5.35E+02	5.36E+02	5.36E+02	5.37E+02	5.37E+02	5.38E+02	5.38E+02	5.39E+02	5.39E+02	5.40E+02	5.40E+02	5.41E+02	5.41E+02	5.42E+02	5.42E+02	5.42E+02	5.42E+02	5.42E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02	5.43E+02];
CO225data = [2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.78E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.79E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.80E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.81E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.82E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.83E+02	2.83E+02	2.83E+02	2.83E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.84E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.85E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.86E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.87E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.88E+02	2.89E+02	2.89E+02	2.89E+02	2.90E+02	2.90E+02	2.91E+02	2.91E+02	2.92E+02	2.92E+02	2.93E+02	2.93E+02	2.93E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.94E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.95E+02	2.96E+02	2.96E+02	2.96E+02	2.96E+02	2.97E+02	2.97E+02	2.98E+02	2.98E+02	2.99E+02	2.99E+02	2.99E+02	3.00E+02	3.00E+02	3.00E+02	3.01E+02	3.01E+02	3.01E+02	3.02E+02	3.02E+02	3.02E+02	3.03E+02	3.03E+02	3.03E+02	3.04E+02	3.04E+02	3.05E+02	3.05E+02	3.05E+02	3.06E+02	3.06E+02	3.07E+02	3.07E+02	3.08E+02	3.08E+02	3.09E+02	3.09E+02	3.09E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.10E+02	3.11E+02	3.11E+02	3.11E+02	3.12E+02	3.12E+02	3.12E+02	3.13E+02	3.14E+02	3.14E+02	3.15E+02	3.16E+02	3.16E+02	3.17E+02	3.18E+02	3.18E+02	3.19E+02	3.20E+02	3.21E+02	3.22E+02	3.23E+02	3.24E+02	3.25E+02	3.26E+02	3.27E+02	3.29E+02	3.30E+02	3.31E+02	3.32E+02	3.33E+02	3.35E+02	3.37E+02	3.38E+02	3.40E+02	3.41E+02	3.42E+02	3.44E+02	3.45E+02	3.47E+02	3.49E+02	3.51E+02	3.52E+02	3.54E+02	3.55E+02	3.56E+02	3.57E+02	3.58E+02	3.60E+02	3.61E+02	3.63E+02	3.65E+02	3.67E+02	3.69E+02	3.70E+02	3.73E+02	3.75E+02	3.77E+02	3.79E+02	3.81E+02	3.83E+02	3.85E+02	3.87E+02	3.89E+02	3.92E+02	3.94E+02	3.96E+02	3.98E+02	4.01E+02	4.03E+02	4.05E+02	4.08E+02	4.10E+02	4.12E+02	4.14E+02	4.17E+02	4.19E+02	4.21E+02	4.23E+02	4.24E+02	4.26E+02	4.28E+02	4.29E+02	4.31E+02	4.32E+02	4.33E+02	4.35E+02	4.36E+02	4.37E+02	4.38E+02	4.38E+02	4.39E+02	4.40E+02	4.40E+02	4.41E+02	4.41E+02	4.41E+02	4.42E+02	4.42E+02	4.42E+02	4.42E+02	4.42E+02	4.43E+02	4.43E+02	4.43E+02	4.43E+02	4.43E+02	4.43E+02	4.43E+02	4.42E+02	4.42E+02	4.42E+02	4.42E+02	4.42E+02	4.41E+02	4.41E+02	4.41E+02	4.40E+02	4.40E+02	4.40E+02	4.39E+02	4.39E+02	4.38E+02	4.37E+02	4.37E+02	4.36E+02	4.36E+02	4.35E+02	4.35E+02	4.34E+02	4.33E+02	4.33E+02	4.32E+02	4.32E+02	4.31E+02	4.31E+02	4.30E+02	4.29E+02	4.29E+02	4.28E+02	4.28E+02	4.27E+02	4.27E+02	4.26E+02	4.25E+02	4.25E+02	4.24E+02	4.24E+02	4.23E+02	4.23E+02	4.22E+02	4.22E+02	4.21E+02	4.21E+02	4.20E+02	4.20E+02	4.19E+02	4.19E+02	4.19E+02	4.18E+02	4.18E+02	4.17E+02	4.17E+02	4.16E+02	4.16E+02	4.15E+02	4.15E+02	4.14E+02	4.14E+02	4.14E+02	4.13E+02	4.13E+02	4.12E+02	4.12E+02	4.11E+02	4.11E+02	4.10E+02	4.10E+02	4.09E+02	4.09E+02	4.08E+02	4.08E+02	4.08E+02	4.07E+02	4.07E+02	4.06E+02	4.06E+02	4.05E+02	4.05E+02	4.04E+02	4.04E+02	4.04E+02	4.03E+02	4.03E+02	4.03E+02	4.02E+02	4.02E+02	4.01E+02	4.01E+02	4.01E+02	4.00E+02	4.00E+02	3.99E+02	3.99E+02	3.99E+02	3.98E+02	3.98E+02	3.98E+02	3.97E+02	3.97E+02	3.97E+02	3.96E+02	3.96E+02	3.96E+02	3.95E+02	3.95E+02	3.95E+02	3.94E+02	3.94E+02	3.94E+02	3.93E+02	3.93E+02	3.93E+02	3.92E+02	3.92E+02	3.92E+02	3.92E+02	3.91E+02	3.91E+02	3.91E+02	3.90E+02	3.90E+02	3.90E+02	3.89E+02	3.89E+02	3.89E+02	3.89E+02	3.88E+02	3.88E+02	3.88E+02	3.87E+02	3.87E+02	3.87E+02	3.86E+02	3.86E+02	3.86E+02	3.86E+02	3.85E+02	3.85E+02	3.85E+02	3.85E+02	3.84E+02	3.84E+02	3.84E+02	3.83E+02	3.83E+02	3.83E+02	3.83E+02	3.82E+02	3.82E+02	3.82E+02	3.82E+02	3.81E+02	3.81E+02	3.81E+02	3.81E+02	3.80E+02	3.80E+02	3.80E+02	3.80E+02	3.79E+02	3.79E+02	3.79E+02	3.79E+02	3.78E+02	3.78E+02	3.78E+02	3.78E+02	3.77E+02	3.77E+02	3.77E+02	3.77E+02	3.76E+02	3.76E+02	3.76E+02	3.76E+02	3.75E+02	3.75E+02	3.75E+02	3.75E+02	3.74E+02	3.74E+02	3.74E+02	3.74E+02	3.74E+02	3.73E+02	3.73E+02	3.73E+02	3.73E+02	3.72E+02	3.72E+02	3.72E+02	3.72E+02	3.71E+02	3.71E+02	3.71E+02	3.71E+02	3.70E+02	3.70E+02	3.70E+02	3.70E+02	3.70E+02	3.69E+02	3.69E+02	3.69E+02	3.69E+02	3.69E+02	3.68E+02	3.68E+02	3.68E+02	3.68E+02	3.67E+02	3.67E+02	3.67E+02	3.67E+02	3.67E+02	3.66E+02	3.66E+02	3.66E+02	3.66E+02	3.65E+02	3.65E+02	3.65E+02	3.65E+02	3.65E+02	3.64E+02	3.64E+02	3.64E+02	3.64E+02	3.64E+02	3.63E+02	3.63E+02	3.63E+02	3.63E+02	3.63E+02	3.62E+02	3.62E+02	3.62E+02	3.62E+02	3.62E+02	3.61E+02	3.61E+02	3.61E+02	3.61E+02	3.60E+02	3.60E+02	3.60E+02	3.60E+02	3.60E+02	3.60E+02	3.59E+02	3.59E+02	3.59E+02	3.59E+02	3.58E+02	3.58E+02	3.58E+02	3.58E+02	3.58E+02	3.58E+02	3.57E+02	3.57E+02	3.57E+02	3.57E+02	3.57E+02	3.56E+02	3.56E+02	3.56E+02	3.56E+02	3.56E+02	3.55E+02	3.55E+02	3.55E+02	3.55E+02	3.55E+02	3.54E+02	3.54E+02	3.54E+02	3.54E+02	3.54E+02	3.54E+02	3.53E+02	3.53E+02	3.53E+02	3.53E+02	3.53E+02	3.52E+02	3.52E+02	3.52E+02	3.52E+02	3.52E+02	3.52E+02	3.51E+02	3.51E+02	3.51E+02	3.51E+02	3.51E+02	3.50E+02	3.50E+02	3.50E+02	3.50E+02	3.50E+02	3.50E+02	3.49E+02	3.49E+02	3.49E+02	3.49E+02	3.49E+02	3.48E+02	3.48E+02	3.48E+02	3.48E+02	3.48E+02	3.48E+02	3.47E+02	3.47E+02	3.47E+02	3.47E+02	3.47E+02	3.47E+02	3.46E+02	3.46E+02	3.46E+02	3.46E+02	3.46E+02	3.46E+02	3.45E+02	3.45E+02	3.45E+02	3.45E+02	3.45E+02	3.44E+02	3.44E+02	3.44E+02	3.44E+02	3.44E+02	3.44E+02	3.44E+02	3.43E+02	3.43E+02	3.43E+02	3.43E+02	3.43E+02	3.42E+02	3.42E+02	3.42E+02	3.42E+02	3.42E+02	3.42E+02	3.42E+02	3.41E+02	3.41E+02	3.41E+02	3.41E+02	3.41E+02	3.41E+02	3.40E+02	3.40E+02	3.40E+02	3.40E+02	3.40E+02	3.40E+02	3.39E+02	3.39E+02	3.39E+02	3.39E+02	3.39E+02	3.39E+02	3.38E+02	3.38E+02	3.38E+02	3.38E+02	3.38E+02	3.38E+02	3.37E+02	3.37E+02	3.37E+02	3.37E+02	3.37E+02	3.37E+02	3.37E+02	3.36E+02	3.36E+02	3.36E+02	3.36E+02	3.36E+02	3.36E+02	3.35E+02	3.35E+02	3.35E+02	3.35E+02	3.35E+02	3.35E+02	3.35E+02	3.34E+02	3.34E+02	3.34E+02	3.34E+02	3.34E+02	3.34E+02	3.33E+02	3.33E+02	3.33E+02	3.33E+02	3.33E+02	3.33E+02	3.33E+02	3.32E+02	3.32E+02	3.32E+02	3.32E+02	3.32E+02	3.32E+02	3.32E+02	3.31E+02	3.31E+02	3.31E+02	3.31E+02	3.31E+02	3.31E+02	3.30E+02	3.30E+02	3.30E+02	3.30E+02	3.30E+02	3.30E+02	3.30E+02	3.29E+02	3.29E+02	3.29E+02	3.29E+02	3.29E+02	3.29E+02	3.29E+02	3.28E+02	3.28E+02	3.28E+02	3.28E+02	3.28E+02	3.28E+02	3.28E+02	3.27E+02	3.27E+02	3.27E+02];

%reshaped cause I like to look at columns as opposed to rows 
CO225data = reshape(CO225data, [length(CO285data) 1]);
CO245data = reshape(CO245data, [length(CO285data) 1]);
CO260data = reshape(CO260data, [length(CO285data) 1]);
CO285data = reshape(CO285data, [length(CO285data) 1]);
%%
carbon = CO285data;
methane = CH485data;

%BAND PROPERTIES
Q6065 = 225.0452;
Q6570 = 204.9277;
Q7075 = 191.2693; 

A6065 = 1.0273e+13; 
A6570 = 8.5146e+12;
A7075 = 6.6906e+12;
Atot = A6065+A6570+A7075;

f6065 = A6065./Atot;
f6570 = A6570./Atot;
f7075 = A7075./Atot;

%instantiate ebm object and prescribe forcings 
e6065 = ebm(9.75, 104, Q6065, 0.6);
e6570 = ebm(9.75, 104, Q6570, 0.6);
e7075 = ebm(9.75, 104, Q7075, 0.6);

%% Euler method
x0=0.9;

%initial forcings: 
e6065.mu = carbon(1);
e6065.nu = methane(1);
e6570.mu = carbon(1);
e6570.nu = methane(1);
e7075.mu = carbon(1);
e7075.nu = methane(1);

%solve and store tau_s
eqtemp(e6065, x0);
eqtemp(e6570, x0);
eqtemp(e7075, x0);

tau6065(1) = e6065.tau_s;
tau6570(1) = e6570.tau_s;
tau7075(1) = e7075.tau_s;

%sum masses so each has same total CO2 and CH4
mass(1,:) = e6065.state(end-1:end); %could've used any ebm

%update states
%dynamics
de6065 = dynamics(e6065);
m6065 = de6065(end-1:end).*[f6065 f6065];
de6570 = dynamics(e6570);
m6570 = de6570(end-1:end).*[f6570 f6570];
de7075 = dynamics(e7075);
m7075 = de7075(end-1:end).*[f7075 f7075];

M = m6065 + m6570 + m7075;

%account for work above 
D6065 = de6065;
D6065(end-1:end) = M;
D6570 = de6570;
D6570(end-1:end) = M;
D7075 = de7075;
D7075(end-1:end) = M;
e6065.state = e6065.state + D6065;
e6570.state = e6570.state + D6570;
e7075.state = e7075.state + D7075;

for i = 1:length(carbon)-1
    %apply forcings 
    e6065.mu = carbon(i+1);
    e6065.nu = methane(i+1);
    e6570.mu = carbon(i+1);
    e6570.nu = methane(i+1);
    e7075.mu = carbon(i+1);
    e7075.nu = methane(i+1);
    
    %solve and store tau_s
    eqtemp(e6065, e6065.tau_s);
    eqtemp(e6570, e6570.tau_s);
    eqtemp(e7075, e7075.tau_s);

    tau6065(end+1) = e6065.tau_s;
    tau6570(end+1) = e6570.tau_s;
    tau7075(end+1) = e7075.tau_s;
    
    %update state
    %dynamics
    de6065 = dynamics(e6065);
    m6065 = de6065(end-1:end).*[f6065 f6065];
    de6570 = dynamics(e6570);
    m6570 = de6570(end-1:end).*[f6570 f6570];
    de7075 = dynamics(e7075);
    m7075 = de7075(end-1:end).*[f7075 f7075];

    M = m6065 + m6570 + m7075;

    %account for work above 
    D6065 = de6065;
    D6065(end-1:end) = M;
    D6570 = de6570;
    D6570(end-1:end) = M;
    D7075 = de7075;
    D7075(end-1:end) = M;
    e6065.state = e6065.state + D6065;
    e6570.state = e6570.state + D6570;
    e7075.state = e7075.state + D7075;
    
    mass(end+1,:) = e6065.state(end-1:end);
    mass(end)
end