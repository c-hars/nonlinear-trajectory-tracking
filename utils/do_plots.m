%% Draw plots (position, attitude, actuator inputs)

function do_plots(t,X,U,r_,qp)

nr = size(r_,2);

t = t(1:nr);
X = X(1:nr,:);

cmap = colororder;

%%

subplot(3,1,1); hold on

plot(t,X(:,1))
plot(t,X(:,2))
plot(t,X(:,3))

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

if size(X,2) == 13
    % Quaternion
    q_ = X(:,7:10);
    eul = quat2eul(q_);
    plot(t,eul(:,3)*180/pi)
    plot(t,eul(:,2)*180/pi)
    plot(t,eul(:,1)*180/pi)
else
    % Assume 12 state and Euler
    plot(t,X(:,7)*180/pi)
    plot(t,X(:,8)*180/pi)
    plot(t,X(:,9)*180/pi)
end

grid on
yline(0,'Color',[1 1 1]*0.25)
ylabel('Euler angles [deg]')
legend('\phi','\theta','\psi')

xlim([0 t(end)])
ylim([-1 1]*60)

%%

subplot(3,1,3)

plot(t(1:end-1), U(1:nr-1,:)) % control input

grid on
yline(0,'Color',[1 1 1]*0.25)
ylabel('\delta U')

xlim([0 t(end)])
ylim([-1 1]*400)

%%

xlabel('Time [s]')

end