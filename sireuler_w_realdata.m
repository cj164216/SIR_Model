
%%
%Initialization
T = 360;
n = 1000;
h = T/n;
E = zeros(3,T/h);

%%
%The basic information about COVID-19
R0 = 2.75 ; % From the published R0 data with the original variant 病毒基本属性
%R0 = lamda/mu, the definition
mu = 1/14; %recovering rate, generally need 14 days to recover
lamda0 = R0*mu; %initial infection rate(speed),
%lamada: the number of a infected one can infect others in a day without protection
%lamda 会变化，因为封城会降低，
%lamda = ncon*P_inf， ncon是接触人数，P_inf是接触后的感染概率
%ncon可变，P_inf不会变(病毒的基本属性)
ncon0 = 5; %initially suppose a infected one contact with 5 people in a day without any protection
P_inf = lamda0/ncon0; %infection probability,0.03926

% adjustable parameter
%因为防疫政策可以变化
ncon = 2; %number of contacted peple within a day with different restrict
lamda = ncon*P_inf;

%%
%relate to the total population of Wuhan
N = 11e6; %population of Wuhan
E(1,1) = 425/N; %infection
E(3,1) = 28/N; %recover
E(2,1) = 1-E(1,1)-E(3,1); %suspection


%% 
%main loop for Euler method
for i = 1:(T/h)-1
    E(1,i+1) = E(1,i) + (lamda*E(1,i)*E(2,i) - mu*E(1,i))*h;
    E(2,i+1) = E(2,i) + (-lamda*E(1,i)*E(2,i)*h);
    E(3,i+1) = E(3,i) + (mu*E(1,i)*h);
end

%%
%convert from ratio to number of people
c = E(3,:)*N;
b = E(2,:)*N;
a = E(1,:)*N;

%%
%plot part
d = zeros(1,n);
for j= 1:n
    d(:,j) = j*h;
end

figure
hold all;
plot(d,a,'blue');
%plot(d,b,'green');
plot(d,c,'red');