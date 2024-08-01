% *************************************************************************
% CMI2NI: Conditional mutual inclusive information(CMI2)-based Network
% Inference method from gene expression data  
% *************************************************************************
% This is matlab code for netwrk inference method CMI2NI. 
% Input: 输入
% 'data' is expression of variable,in which row is varible and column is the sample;
% 'data'的行是基因 列是样本
% 'lamda' is the parameter decide the dependence;
% 'lamda'是决定独立性的参数
% 'order0' is the parameter to end the program when order=order0;
% 'order0'是结束程序的参数,当order=order0时 结束程序
% If nargin==2,the algorithm will be terminated untill there is no change 
% in network toplogy.
% 如果nargin==2,直到网络拓扑没有变化 程序才会结束
% Output: 输出
% 'G' is the 0-1 network or graph after pc algorithm;
% 'G'是pc算法后的0-1网络或图形
% 'Gval' is the network with strenthness of dependence;
% 'Gval'是具有依赖性强度的网络
% 'order' is the order of the pc algorithm, here is equal to order0;
% 'order'是pc算法的阶数,这里等于order0
% Example:
% 
% Author: Xiujun Zhang.
% Version: Sept.2014.

function [G,Gval,order]=CMI2NI(data,lamda,order0)                                                 
n_gene=size(data,1);% n_gene是基因的数量
G=ones(n_gene,n_gene);% G是元素全为1的矩阵
G=tril(G,-1)';% tril(G,-1)提取主对角线下方的元素,然后转置
G=G+G';%G 主对角为0,其余元素为1的矩阵
Gval=G;%Gval 主对角线的元素为0,其余元素为1的矩阵
order=-1;t=0;
while t==0
     order=order+1;%order=-1+1=0
     if nargin==3
       if order>order0
           order=order-1; % 
           return
       end
     end
    [G,Gval,t]=edgereduce(G,Gval,order,data,t,lamda);%产生权重矩阵
 
     if t==0
          disp('No edge is reduce! Algorithm  finished!');
          break;
     else 
          t=0;
     end
end
   order=order-1; % The value of order is the last order of the algorithm 
end

%% edgereduce  
function [G,Gval,t]=edgereduce(G,Gval,order,data,t,lamda)
G0=G;% G0 主对角线元素为0 其余元素为1的矩阵
%[nrow,ncol]=find(G~=0);
if order==0
    for i=1:size(G,1)%基因循环
        for j=1:size(G,1)%基因循环
            if G(i,j)~=0%如果G(i,j)不等于0
                cmiv=cmi(data(i,:),data(j,:));%计算第i个基因和第j个基因之间的条件互信息
                Gval(i,j)=cmiv;  Gval(j,i)=cmiv;%Gval矩阵元素为CMI值
                if cmiv<lamda
                    G(i,j)=0;G(j,i)=0;%如果互信息小于阈值,矩阵G的对应元素为0
                end
            end
        end
    end
          t=t+1;
else
  for i=1:size(G,1)%基因循环
      for j=1:size(G,1)%基因循环
          if G(i,j)~=0%当G(i,j)不等于0时,即互信息大于阈值
              adj=[] ;
              for k=1:size(G,1)
                  if G(i,k)~=0 && G(j,k)~=0
                      adj=[adj,k];%记录下与i和j都连接的基因
                  end
              end
              if size(adj,2)>=order%size(adj,2)即文中的L,此处判断共同连接基因数目是否大于阶数
                   combntnslist=combntns(adj,order);%combntns(x,m)列举从n个元素中取出m个元素的组合,x是含有n个元素的向量
                   combntnsrow=size(combntnslist,1);   
                   cmiv=0;
                   v1=data(i,:);v2=data(j,:);
                   for k=1:combntnsrow   
                     vcs=data(combntnslist(k,:),:);   
                     a=MI2(v1,v2,vcs) ;%计算CMI2
                     cmiv=max(cmiv,a);%取最大的那个CMI2
                   end
                   Gval(i,j)=cmiv; Gval(j,i)=cmiv;%Gval矩阵的值更新为CMI2的值
                   if cmiv<lamda
                         G(i,j)=0; G(j,i)=0;%CMI2和阈值lamda比较
                   end              
                   t=t+1; 
              end
          end
                      
      end
  end 
%   if t~=0 && ((size(G0,1)*size(G0,2))-sum(sum(G==G0))<20)
%   if t~=0 && ((size(G0,1)*size(G0,2))-sum(sum(G==G0))==0)
%       t = 0; 
%   end               
end
end

%% compute conditional mutual information of x and y 
function cmiv=cmi(v1,v2,vcs)
 if  nargin==2
        c1=det(cov(v1));%cov求协方差,det求行列式
        c2=det(cov(v2));
        c3=det(cov(v1,v2));
        cmiv=0.5*log(c1*c2/c3);%文中等式2 计算互信息 
     elseif  nargin==3
        c1=det(cov([v1;vcs]'));
        c2=det(cov([v2;vcs]'));
        c3=det(cov(vcs'));
        c4=det(cov([v1;v2;vcs]'));
        cmiv=0.5*log((c1*c2)/(c3*c4));%文中等式4 计算条件互信息       
 end
    % cmiv=abs(cmiv);
     if  cmiv==inf 
            cmiv=1.0e+010;
     end
end

% Conditional mutul inclusive information (CMI2)计算CMI2
function r_dmi = MI2(x,y,z)
r_dmi = (cas(x,y,z) + cas(y,x,z))/2;%文中等式7 等式10
end

% x and y are 1*m dimensional vector; z is n1*m dimensional.
function CS = cas(x,y,z)
% x=rand(10,1)';y=rand(10,1)';z=rand(10,2)';
% x,y,z are row vectors;
n1 = size(z,1);
n = n1 +2;

Cov = cov(x);%cov求协方差
Covm = cov([x;y;z]');
Covm1 = cov([x;z]');

InvCov = inv(Cov);%inv是求矩阵的逆
InvCovm = inv(Covm);
InvCovm1 = inv(Covm1);

C11 = InvCovm1(1,1);%Cxx
C12 = 0;%Cxy
C13 = InvCovm1(1,2:1+n1);%Cxz
C23 = InvCovm(2,3:2+n1)-InvCovm(1,2) * (1/(InvCovm(1,1)-InvCovm1(1,1)+InvCov(1,1))) * (InvCovm(1,3:2+n1) - InvCovm1(1,2:1+n1)) ;%Cyz
C22 = InvCovm(2,2)- InvCovm(1,2)^2 * (1/(InvCovm(1,1)-InvCovm1(1,1)+InvCov(1,1)));%Cyy
C33 = InvCovm(3:2+n1,3:2+n1)- (1/(InvCovm(1,1)-InvCovm1(1,1)+InvCov(1,1))) * ((InvCovm(1,3:2+n1)-InvCovm1(1,2:1+n1))'*(InvCovm(1,3:2+n1)-InvCovm1(1,2:1+n1)));%Czz
InvC = [[C11,C12,C13];[C12,C22,C23];[[C13',C23'],C33]];%C
% C = inv(InvC);  

C0 = Cov(1,1) * (InvCovm(1,1) - InvCovm1(1,1) + InvCov(1,1));%C0
CS = 0.5 * (trace(InvC*Covm)+log(C0)-n) ;%等式10的1/2
 
end