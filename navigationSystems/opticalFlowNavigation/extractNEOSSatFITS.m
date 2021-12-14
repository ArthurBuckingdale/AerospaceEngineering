%the purpose of this script is to extract data from the fits header and
%transform it into data to be used for Kalman Filter Testing.

% close all
% clear all

% read the images
pathToImages = '/Users/jun/Documents/NEOSSat Images/neossat5 - nav systems proj/images/';
imds = imageDatastore(pathToImages,'ReadFcn',@fitsread);

% %quick slideshow of all the images
% for i = 1:length(imds.Files)
%     im = readimage(imds,i);
%     imagesc(im)
%     pause(0.1)
% end

% locate the needed info in the fits header
for i = 1:length(imds.Files)
    disp(i)
    % define the keys we need
    rightAscensionPointingKey = 'OBJCTRA';
    declinationKey = 'OBJCTDEC';
    timeStamp = 'A_EXP_S';
    stateVectorKey1 = 'JPOS1_1';
    stateVectorKey2 = 'JPOS1_2';
    stateVectorKey3 = 'JPOS1_3';
    stateVectorKey4 = 'JVEL1_1';
    stateVectorKey5 = 'JVEL1_2';
    stateVectorKey6 = 'JVEL1_3';
    rollKey = 'CMDROL';
    quat0 = 'CMDQ0'; 
    quat1 = 'CMDQ1';
    quat2 = 'CMDQ2';
    quat3 = 'CMDQ3';
   
    %open the fits header.
    info = fitsinfo(imds.Files{i}).PrimaryData.Keywords;
   
    % find each one of the target values
    [idx,~] = find(strcmp(info,rightAscensionPointingKey));
    RA(i) = info(idx,2);
    
    [idx2,~] = find(strcmp(info,declinationKey));
    Dec(i) = info(idx2,2);
   
    [idx3,~] = find(strcmp(info,stateVectorKey1));
    stateVect(i,1) = str2num((info{idx3,2}));
   
    [idx4,~] = find(strcmp(info,stateVectorKey2));
    stateVect(i,2) = str2num((info{idx4,2}));
   
    [idx5,~] = find(strcmp(info,stateVectorKey3));
    stateVect(i,3) = str2num((info{idx5,2}));
   
    [idx6,~] = find(strcmp(info,stateVectorKey4));
    stateVect(i,4) = str2num((info{idx6,2}));
   
    [idx7,~] = find(strcmp(info,stateVectorKey5));
    stateVect(i,5) = str2num((info{idx7,2}));
   
    [idx8,~] = find(strcmp(info,stateVectorKey6));
    stateVect(i,6) = str2num((info{idx8,2}));
   
    [idx9,~] = find(strcmp(info,timeStamp));
    timeStamps(i) = info(idx9,2);
    
    [idx10,~] = find(strcmp(info,rollKey));
    spacecraftroll(i) = info(idx10,2);
    
    [idx11,~] = find(strcmp(info,quat0));
    cmdq0(i) = info(idx11,2);
    
    [idx12,~] = find(strcmp(info,quat1));
    cmdq1(i) = info(idx12,2);
    
    [idx13,~] = find(strcmp(info,quat2));
    cmdq2(i) = info(idx13,2);
    
    [idx14,~] = find(strcmp(info,quat3));
    cmdq3(i) = info(idx14,2);
end

dataTable = table(RA', Dec', stateVect(:,1),stateVect(:,2),stateVect(:,3),...
        stateVect(:,4),stateVect(:,5),stateVect(:,6), timeStamps',spacecraftroll',cmdq0',cmdq1',cmdq2',cmdq3',...
        'VariableNames',{'Pointing RA(Hrs Min Sec)','Pointing Dec(Deg Min Sec)','pos_x(Km)',...
        'pos_y(Km)','pos_z(Km)','vel_x(Km/s)','vel_y(Km/s)','vel_z(Km/s)','DateTime',...
        'Spacecraft Roll(Deg)','quat0','quat1','quat2','quat3'});
    
    
writetable(dataTable,'systemesDeNavigation.csv')

for i = 1:length(timeStamps)
    time(i) = str2num(timeStamps{i}(end-11:end-10)) + str2num(timeStamps{i}(end-8:end-7))/60 + str2num(timeStamps{i}(end-5:end))/3600;
    svNorm(i) = norm(stateVect(i,1:3));
    svVelNorm(i) = norm(stateVect(i,4:6));
end


for i = 1:length(RA)
    decDeg(i) =  str2num(Dec{i}(1:2)) + str2num(Dec{i}(4:5))/60 + str2num(Dec{i}(7:9))/3600;
    raDeg(i) = (360 * str2num(RA{i}(1:2)) / 24) + str2num(RA{i}(4:5))*(0.25) + str2num(RA{i}(7:9))*(0.25/60);
    roll(i) = spacecraftroll{i};
    %time(i) = timeStamps{i};
end


figure
hold on
plot(time,stateVect(:,1),'r*')
plot(time,stateVect(:,2),'b*')
plot(time,stateVect(:,3),'g*')
ylabel('Position (Km)')
yyaxis right
plot(time,svVelNorm,'*-')
ylabel('Pos Magnitude')
hold off
legend('X-Pos','Y-Pos','Z-Pos')
xlabel('Observation (Increasing Time)')

figure
hold on
plot(time,stateVect(:,1),'r*')
plot(time,stateVect(:,2),'b*')
plot(time,stateVect(:,3),'g*')
ylabel('Position (Km)')
yyaxis right
plot(time,decDeg,'k+')
ylabel('Declination')
hold off
legend('X-Pos','Y-Pos','Z-Pos')
xlabel('Observation (Increasing Time)')


figure
hold on
plot(time,stateVect(:,4),'r*')
plot(time,stateVect(:,5),'b*')
plot(time,stateVect(:,6),'g*')
hold off
legend('X-Vel','Y-Vel','Z-Vel')
xlabel('Observation (Increasing Time)')
ylabel('Speed (Km/s)')



figure
subplot(1,3,1)
plot(time,decDeg,'k+')
title('Declination')
xlabel('Observation (Increasing Time)')
ylabel('Degrees')
subplot(1,3,2)
plot(time,raDeg,'ko')
title('Right Ascension')
xlabel('Observation (Increasing Time)')
ylabel('Degrees')
subplot(1,3,3)
plot(time,roll,'kx')
title('Roll')
xlabel('Observation (Increasing Time)')
ylabel('Degrees')

% figure
% hold on
% plot(str2num(val),'rx-')
% plot(str2num(val2),'gx-')
% plot(str2num(val3),'bx-')
% hold off

figure
hold on
plot(time,[cmdq0{:}]','*')
plot(time,[cmdq1{:}]','*')
plot(time,[cmdq2{:}]','*')
plot(time,[cmdq3{:}]','*')
hold off
xlabel('Time(hours)')

for i = 1:length(cmdq0)
    quaternion = [cmdq0{i},cmdq1{i},cmdq1{i},cmdq2{i}];
    betterQuaternion(i,:) = [cmdq3{i},cmdq0{i},cmdq1{i},cmdq2{i}];
    % we now want the attitude angles
    quatPsi(i) = atan2d(2*(quaternion(2)*quaternion(3)+quaternion(1)*quaternion(4)),1-2*(quaternion(3).^2+quaternion(4).^2));
    quatTheta(i) = -asind(2*(quaternion(2)*quaternion(4) - quaternion(1)*quaternion(3)));
    quatPhi(i) = atan2d(2*(quaternion(1)*quaternion(2)+quaternion(3)*quaternion(4)),1-2*(quaternion(2).^2+quaternion(3).^2));
end

figure
plot(betterQuaternion)

figure
hold on
plot(time,quatTheta,'*')
plot(time,quatPsi,'*')
plot(time,quatPhi,'*')
hold off


