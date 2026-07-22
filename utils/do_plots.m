%% Draw plots (position, attitude, actuator inputs)

function do_plots(t, X, U, r_, qp, AttRep)

nr = size(r_,2);
t = t(1:nr);
X = X(1:nr,:);
Xe = to_euler_state(X, AttRep);
cmap = colororder;

%%

subplot(3,1,1); hold on

plot(t,Xe(:,1))
plot(t,Xe(:,2))
plot(t,Xe(:,3))

plot(t,r_(1,:),'Color',cmap(1,:),'LineStyle','--')
plot(t,r_(2,:),'Color',cmap(2,:),'LineStyle','--')
plot(t,r_(3,:),'Color',cmap(3,:),'LineStyle','--')

grid on
yline(0,'Color',[1 1 1]*0.25)
ylabel('Position [m]')
legend('x','y','z')

xlim([0 t(end)])
ylim([-1 1]*1.1*max(abs(r_),[],'all'))

%%

subplot(3,1,2); hold on

sas = @(x) (sum(x(x>0)) - sum(abs(x(x<0)))) / (sum(x(x>0)) + sum(abs(x(x<0))));
plot(t, Xe(:,7)*180/pi)
plot(t, Xe(:,8)*180/pi)
plot(t, Xe(:,9)*180/pi)

grid on
yline(0,'Color',[1 1 1]*0.25,'HandleVisibility','off')
ylabel('Euler angles [deg]')
legend('\phi','\theta','\psi')

xlim([0 t(end)])
ylim([-1 1]*60)

%%

subplot(3,1,3); hold on

% Actuator limits + nominal
yline(qp.max_du,'Color',[1 1 1]*0.5)
yline(0,'Color',[1 1 1]*0.67)
yline(-qp.nominal_omegas,'Color',[1 1 1]*0.5)

plot(t(1:end-1), U(1:nr-1,:)) % control input

grid on
ylabel('\delta U')

xlim([0 t(end)])
ylim([-1 1]*400)


%%

xlabel('Time [s]')

end