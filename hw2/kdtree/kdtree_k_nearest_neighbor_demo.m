clc;

%% compile
% mex kdtree_build.cpp
mex kdtree_k_nearest_neighbors.cpp
% disp('compiled.');

%% create data and execute query
%rand('seed',1);
%p = rand( 30, 2 ); % input data
train = load( 'datas/hw3_train.dat' );
p = train( :, 1:2 );
q = [.5,.5]'; % query data
%k = 1; % query how many neighbor
tree = kdtree_build( p );
idxs = kdtree_k_nearest_neighbors(tree,q,k);

%% visualize
close all;
xlim( [0 1] );
ylim( [0 1] );
hold on; %axis equal; axis off;
plot( p(:,1), p(:,2), '.b');
plot(q(1), q(2),'xk');
plot(p(idxs,1), p(idxs,2),'or');
legend('database', 'query', 'query results');
xlabel('x coordinate');
ylabel('y coordinate');
title('3.3(1)');
saveas(gcf, 'outputs/3.3(1).png', 'png');

return;
%% text lables
for i=1:size(q,1)
    text(q(i,1), q(i,2),sprintf(' %d',i));
    text(p(idxs(i),1), p(idxs(i),2),sprintf(' %d',i));
end
