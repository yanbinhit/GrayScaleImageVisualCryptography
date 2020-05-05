prevDir = pwd;
[dir, dummy, dummy2] = fileparts(mfilename('fullpath'));
addpath(genpath(dir),'-begin');
cd(dir);