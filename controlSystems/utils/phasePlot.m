function phasePlot(initialCond, timeSpan, funHandle)
init = initialCond; % Initial Conditions vector
initt = timeSpan(1); finalt = timeSpan(2); % Timespan
butt=1;
figure(1); hold on; grid on;
for i = 1:size(initialCond,1)
    disp(i)
    [T, Y]=ode45(funHandle,[initt finalt],initialCond(i,:));
    %         figure(2);
    %         plot(T,Y); xlabel('time');title('Time Response');
    plot(Y(:,1),Y(:,2),'b');xlabel('x');ylabel('xdot');title('Phase Plane');
    ylim([-2,2])
    xlim([-2,2])
    x1 = gradient(Y(:,1));
    x2 = gradient(Y(:,2));
    quiver(Y(:,1),Y(:,2),x1,x2,'b');
end
hold off
end