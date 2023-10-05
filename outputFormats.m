clearvars
close all

%For exercise-0
fprintf('%8s%9s%14s\n','NumElem','h','Error')  %before the main loop on div. number
err = norm(u-U(nodes(:,1)), inf);              %inside the main loop on div. number
fprintf('%8d%14.6e%14.6e\n',div,h,err)         %inside the main loop on div. number

%For exercise-1 and exercise-2
fprintf('%8s%9s%14s%14s\n','NumNod','x','U','Q')
solution = [(1:numNodes)',nodes,u,Q];
fprintf('%8d%14.6e%14.6e%14.6e\n',solution')