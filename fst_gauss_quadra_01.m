clear all
clc
% element=[1 2;2 3;2 4];
% numelem=size(element,1);
% numnode=4;
% K=zeros(numnode,numnode);
% 
% for e=1:numelem ;
% index=element(e,:) ;
% k=[1 -1;-1 1];
% K(index,index)=K(index,index)+k
% end
% clc
% Le = 1/4;
% syms x
% N = [  1-(x/Le)   4*x  ]
% U = [  0.4216  ;  0.4463  ]
% 
% 
% N(1) = N(1) + 1
% N(2) = N(2) - 1
% 
% N(1) = N(1) + 1 
% N(2) = N(2) - 1
% 
% N(1) = N(1) + 1 
% N(2) = N(2) - 1

% A = [x ,2*x,4 ,5]
% diff(A)
% 
% f = input('Enter your function: ')
% a = input('lower limit: ')
% b = input('upper limit: ')
% g = input('gauss point: ')
% 
% if g==2 
%     w1 = 1; 
%     w2 = 1;
%     x1 = 1/sqrt(3); 
%     x2 = -1/sqrt(3);
%     F1 = ((w1*f(x1)) + (w2*f(x2)))
% else
%     w1 = 5/9; w2 = 8/9; w3 = 5/9
%     x1 = sqrt(3/5); x2 = 0 ; x3 = -sqrt(3/5)
%     F2 = ((w1*f(x1)) + (w2*f(x2)) + (w3*f(x3)))
% 
% end
