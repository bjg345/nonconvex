id = str2num(getenv('SGE_TASK_ID'));

x = readmatrix('x.csv');
x(1:100)

y = readmatrix('y.csv');
y(1:100, 1:2)
