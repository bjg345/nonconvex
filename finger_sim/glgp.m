id = str2num(getenv('SGE_TASK_ID'));

grid_train = readmatrix('./grid.train.csv');
grid_train = readmatrix('./grid.train.csv');
ground_train = readmatrix('/ground.train.csv');
ground_test = readmatrix('./ground.test.csv');
grid_test = readmatrix('./grid.test.csv');

if id <= 500
    stdev = .1;
elseif id <= 1000
    stdev = .25;
else stdev = 1;
end

n = 250;

rng(id);

ind_train = randsample(10000, .8*n);
ind_test = randsample(10000, .2*n);

noise_train = normrnd(0, stdev, .8*n, 1);
noise_test = normrnd(0, stdev, .2*n, 1);

X1 = grid_train(ind_train, :);
X = grid_test(ind_test, :);
y = ground_train(ind_train) + noise_train;
f = ground_test(ind_test);

param  = fit(X1, X, y);
vals_pred = CovMatrix_predict(X1, X, y, param(1), param(2), param(3), param(4));

err = vals_pred-f;
vals_test = f+noise_test;

save(strcat('glgp/err_', num2str(n), '_', num2str(stdev), '_', num2str(id), '.mat'), 'err')
save(strcat('glgp/fit_', num2str(n), '_', num2str(stdev), '_', num2str(id), '.mat'), 'param')
save(strcat('glgp/pred_', num2str(n), '_', num2str(stdev), '_', num2str(id), '.mat'), 'vals_pred')
save(strcat('glgp/vals_', num2str(n), '_', num2str(stdev), '_', num2str(id), '.mat'), 'vals_test')


