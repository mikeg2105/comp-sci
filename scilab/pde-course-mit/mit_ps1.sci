
function probs1

x=0:0.1:2;
fx=1/2*(1-x);
figure(1),plot(x,fx,'k',x-2,fx,'k',x-4,fx,'k',x+2,fx,'k',[-3,3],[0,0],'k--',...
    -2:2:2,zeros(size(-2:2:2)),'k.',-2:2:2,0.5*ones(size(-2:2:2)),'ko',-2:2:2,-0.5*ones(size(-2:2:2)),'ko');
xlabel('x','fontname','times','fontsize',14);
set(gca,'fontname','times','fontsize',14,'xlim',[-3,3],'ylim',[-0.5,0.5]);
set(gca,'plotboxaspectratio',[1,0.3,1]);
%saveas(1,'sine.eps');

fx=abs(1/2*(1-x));
figure(2),plot(x,fx,'k',x-2,fx,'k',x-4,fx,'k',x+2,fx,'k');
xlabel('x','fontname','times','fontsize',14);
set(gca,'fontname','times','fontsize',14,'xlim',[-3,3],'ylim',[0,0.5]);
set(gca,'plotboxaspectratio',[1,0.3,1]);
//saveas(2,'cosine.eps');


x=0:0.01:1; t0=0;
u0=uxt(x,0);
u02=uxt(x,0.2);
u05=uxt(x,0.5);
u1=uxt(x,1);

figure(3),plot(x,u0,'k--',x,u02,'k',x,u05,'k-.',x,u1,'k');
xlabel('x','fontname','times','fontsize',14);
ylabel('u(x,t_0)','fontname','times','fontsize',14);
endfunction
