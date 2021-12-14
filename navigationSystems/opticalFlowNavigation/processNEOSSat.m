% the purpose of this script is to process the NEOSSat images for the
% navigation project. There are a two different things that will be
% extracted here:

close all
clear all

% read the directory and store them
pathToImages = '/Users/jun/Documents/NEOSSat Images/neossat5 - nav systems proj/images';
imds = imageDatastore(pathToImages,'ReadFcn',@fitsread);

% read the dark frame for subtracting
darkFrame = fitsread('/Users/jun/Documents/NEOSSat Images/neossat5 - nav systems proj/aaa-Dark-Frame/NEOS_SCI_2017284084010.fits');
darkFrame = rescale(darkFrame);

%preparing the meta data from the fits file.
for i = 1:length(imds.Files)
    if mod(i,10) == 0
        fprintf('Fits Reading Iteration %d \n',i)
    end
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



% begin the for loop for processing
% this for loop
for i = 1:length(imds.Files)
    if mod(i,10) == 0
        fprintf('Iteration: %d \n',i)
        fprintf('File Name is %s \n',imds.Files{i}(end-26:end))
    end
    % image pre processing
    image = readimage(imds,i);
    image = rescale(image, 0, 1);
    rawIm = image;
    image = image - darkFrame;
    image = rescale(image,0,1);
    image = medfilt2(image, [3,3]);
    binIm = imbinarize(image,'adaptive');
    
%     if mod(i,20) == 0
%         figure
%         subplot(1,2,1)
%         imagesc(binImWrong)
%         title('Binarized with 3x3 medfilt')
%         subplot(1,2,2)
%         imagesc(binIm)
%         title('Binarized with 2x2 medfilt')
%         axis equal
%         linkaxes
%        
%     end
    
    % recovering the stars from the frame
    qq = regionprops(binIm,image,{'Centroid', 'BoundingBox','MaxIntensity'});
    dataVault(i).stars = [cat(1,qq.Centroid),cat(1,qq.BoundingBox),cat(1,qq.MaxIntensity)];
    if mod(i,40) == 0 
        disp(dataVault(i).stars)
        figure
        hold on
        imagesc(binIm)
        plot(dataVault(i).stars(:,1),dataVault(i).stars(:,2),'r*')
        hold off
        title('Binarized Image With Centroids Displayed')
        
    end
    
end


% we now have the dataVault.stars which contains the position of all stars
% present on the image frame. We now need to create a pattern and match the
% stars together with a subsequent frame. To do this, we'll grab our star
% pattern matching algorithm.


%%%%%%%%%%%%%%%%% I cannot send the code for star pattern matching as it
%%%%%%%%%%%%%%%%% belongs to the company I work for. I have included it
%%%%%%%%%%%%%%%%% here as a .P file. 


try
for i = 1:length(dataVault)-1
    if time(i+1)-time(i) > (1/10)
        fprintf('Images %d and %d seperated by more than 1/20 hours \n',i,i+1)
        continue
    end
    coordsImageOne = dataVault(i).stars(:,1:2);
    coordsImageTwo = dataVault(i+1).stars(:,1:7);
    %coordsImageOne = coordsImageOne(coordsImageOne(:,7)>0.6,1:2);
    coordsImageTwo = coordsImageTwo(coordsImageTwo(:,7)>0.55,1:2);
    fprintf('Number of 1st im stars: %d, Number of 2nd Im stars is %d \n',length(coordsImageOne),length(coordsImageTwo))
    index_star(i).idx = star_pattern_matching_proj(coordsImageOne, coordsImageTwo);
    if length(index_star(i).idx) < 4
        continue
    end
    deltas = [];
    for j = 1:length(index_star(i).idx)
        fprintf('Star number %d of iteration %d \n',j,i)
        matchedStars(j,:) = [coordsImageTwo(j,:),coordsImageOne(index_star(i).idx(j),:)];
        deltas(j,:) = norm(matchedStars(j,1:2) - matchedStars(j,3:4));
    end
    tmp = deltas > 200;
    matchToleranceStars(i).matches = matchedStars(~tmp,:);
        if mod(i,25) == 0
            figure
            hold on
            plot(matchToleranceStars(i).matches(:,1),matchToleranceStars(i).matches(:,2),'*')
            plot(matchToleranceStars(i).matches(:,3),matchToleranceStars(i).matches(:,4),'*')
            hold off
            title('Matched Stars In Between Two Sets of Images')
            legend('First Image','Second Image')
            figure
            plot(deltas)
            xlabel('Matched Star')
            ylabel('Norm Pixel Difference')
            title('Residuals of Matched Stars Between Two Images')
        end
    fprintf('Number of matched stars is %d \n',length(matchToleranceStars(i).matches))
end
catch
    
end


% we now need to estimate the pixel delta between each image to obtain the
% velocity estimate. See the algorithm in the LaTeX document to see how all
% the following nonsense works. There are some interesting items here, we
% can used the theta_recovered in order to estimate the roll into the next
% frame. We'll need to keep track of the roll.
clear translation_X_recovered
clear translation_Y_recovered
clear theta_recovered
clear normTrans
clear stateTime
clear newTime
clear deltaVect
n=1;
for i = 1:length(matchToleranceStars)
    if length(matchToleranceStars(i).matches) > 5
        trans = fitgeotrans(matchToleranceStars(i).matches(:,1:2),matchToleranceStars(i).matches(:,3:4),'nonreflectivesimilarity');
        sss = trans.T(2, 1);
        scc = trans.T(1, 1);
        scale_recovered = sqrt(sss*sss+scc*scc);
        translation_X_recovered(n) = trans.T(3, 1);
        translation_Y_recovered(n) = trans.T(3, 2);
        normTrans(n) = norm([translation_X_recovered(n),translation_Y_recovered(n)]);
        theta_recovered(n) = atan2(sss, scc) * 180 / pi; %[deg]
        newTime(n) = time(i);
        stateTime(n,:) = stateVect(i,1:6);
        deltaVect(n,:) = mean(matchToleranceStars(i).matches(:,1:2) - matchToleranceStars(i).matches(:,3:4));
        deltaVectIndividual(n).val = matchToleranceStars(i).matches(:,1:2) - matchToleranceStars(i).matches(:,3:4);
        matchedCoords(n).val = matchToleranceStars(i).matches(:,1:2);
        fprintf('The X and Y translations recovered are [%d,%d] \n',translation_X_recovered(n),translation_Y_recovered(n))
        n=n+1;
    end
end


figure
hold on
plot(newTime, deltaVect(:,1),'*')
plot(newTime, deltaVect(:,2),'*')
hold off
title('Mean Pixel Deltas in each dimension for each image passing filters')
xlabel('Time (h)')
ylabel('Pixel Delta (pixels)')

fprintf('Number of Images giving a velocity measure is %d \n',n)
figure
hold on
plot(newTime, stateTime(:,4),'o')
plot(newTime, stateTime(:,5),'o')
plot(newTime, stateTime(:,6),'o')
ylabel('Velocity Components (km/s)')
legend('X','Y','Z')
yyaxis right
plot(newTime,normTrans,'*')
ylabel('Pixel Shift ')
hold off
title('State Vector and Image Recovered Pixel Velocity Magnitude')

disp('Velocities Measured Code Complete')

% now we need to take the pixel velocities and estimate our linear
% velocities from them. For the sake of not losing my mind, and because i
% cannot un-fuck the rotation matrices. Couple of hardware parameters
% needed here as well. See the NEOSSat specifications for more details. We
% also need to assume a couple of things about the stars present in the
% image. We're going to assume their fixed in the background and that their
% distance is >>> than focal length. We're also going to assume that all
% stars are 10 light years away from us(we are only measuring bright stars,
% so this is not that rediculous of an assumption.

% i'm going to place the coordinates into the frame of

[ra,dec,range]=cart2sph(stateTime(1,1),stateTime(1,2),stateTime(1,3));
[velRa, velDec, velRange] = cart2sph(stateTime(1,4),stateTime(1,5),stateTime(1,6));

radec = [ra,dec];
velradec = [velRa,velDec];

test = dot(radec/norm(radec),velradec/norm(velradec));

fprintf('ra: %d, dec: %d \n',ra, dec)
fprintf('velRa: %d, velDec %d \n', velRa,velDec)
fprintf('pointRA: %d, pointDec: %d \n',raDeg(1) * pi/180,decDeg(1)* pi/180)

q=1;
clear x finalTime groundTruthVelocity groundTruthAng resnorm xDataFilt yDataFilt

for i = 1:size(deltaVect,1)-1
    if abs(newTime(i) - newTime(i+1)) > 1/20
        fprintf('Images [%d,%d] Too Far Apart to Share Stars \n',i,i+1)
        continue
    end
    
    clear xDataFilt
    clear yDataFilt
    
    xData =  512 - matchedCoords(i).val;
    yData = deltaVectIndividual(i).val/((newTime(i) - newTime(i+1))*3600);
    %filter bad data from yData
    yDataMean = median(yData,1);
    qq=1;
    for ff = 1:size(yData,1)
        if abs((yData(ff,1)-yDataMean(1))/yData(ff,1)) < 0.1 & abs((yData(ff,2)-yDataMean(2))/yData(ff,2)) < 0.1
            yDataFilt(qq,:) = yData(ff,:);
            xDataFilt(qq,:) = xData(ff,:);
            qq=qq+1;
        end
    end
    xInitial = [stateTime(i,4:6),0.000005,0.000005,0.000005];
    lb(1:3) = -norm(stateTime(i,4:6));
    ub(1:3) = norm(stateTime(i,4:6));
    options = optimoptions('lsqcurvefit','Display','off',...
        'MaxIterations',1000,'FunctionTolerance',1e-18,'StepTolerance',1e-18);
    [x(q,:),resnorm(q),residual,exitflag,output,lambda,jacobian]  = lsqcurvefit(@computePixelVelFromCameraVel,xInitial,xDataFilt,yDataFilt,[],[],options);
    tangentialVelocity(q,:) = cross(x(q,4:6),stateTime(i,1:3));
    finalTime(q) = newTime(i);
    groundTruthVelocity(q,:) = stateTime(i,4:6);
    groundTruthAng(q,:) = cross(stateTime(i,1:3), stateTime(i,4:6))./(norm(stateTime(i,1:3)).^2);
    q=q+1;
end

figure
plot(resnorm)
title('Residual Normals for Each of the Curve Fitting Operations')
xlabel('Velocity Measure Number')
ylabel('Residual Normal Value')


% figure
% hold on
% plot(timeStamp,spacecraftVelocity(:,1),'*')
% plot(timeStamp,spacecraftVelocity(:,2),'*')
% hold off

% figure
% subplot(2,2,1)
% hold on
% plot(finalTime,x(:,1),'*')
% plot(finalTime,groundTruthVelocity(:,1),'o')
% plot(finalTime, tangentialVelocity(:,1),'+')
% hold off
% subplot(2,2,2)
% hold on
% plot(finalTime,x(:,2),'*')
% plot(finalTime,groundTruthVelocity(:,2),'o')
% plot(finalTime, tangentialVelocity(:,2),'+')
% hold off
% subplot(2,2,3)
% hold on
% plot(finalTime,x(:,3),'*')
% plot(finalTime,groundTruthVelocity(:,3),'o')
% plot(finalTime, tangentialVelocity(:,3),'+')
% hold off
% subplot(2,2,4)
% hold on
% plot(finalTime,vecnorm(x(:,1:3),2,2),'*')
% plot(finalTime,vecnorm(groundTruthVelocity(:,1:3),2,2),'o')
% plot(finalTime,resnorm,'+')
% hold off
% xlabel('Time (hours)')
% ylabel('Velocity m/s')
% legend({'$\hat{x}$','$\hat{y}$','$\hat{z}$','x','y','z'},'Interpreter','latex')


figure
subplot(2,2,1)
hold on
%plot(finalTime,x(:,4),'*')
plot(finalTime, tangentialVelocity(:,1),'+')
ylabel('Computed Velocity (m/s)')
yyaxis right
plot(finalTime,groundTruthVelocity(:,1),'o')
hold off
title({'$v_x$ as a function of time'},'Interpreter','latex')
ylabel('Velocity (m/s)')
xlabel('Time (h)')
legend('Computed Tangential Velocity','Ground Truth Velocity')
subplot(2,2,2)
hold on
%plot(finalTime,x(:,5),'*')
plot(finalTime, tangentialVelocity(:,2),'+')
ylabel('Computed Velocity (m/s)')
yyaxis right
plot(finalTime,groundTruthVelocity(:,2),'o')
hold off
title({'$v_y$ as a function of time'},'Interpreter','latex')
ylabel('Velocity (m/s)')
xlabel('Time (h)')
legend('Computed Tangential Velocity','Ground Truth Velocity')
subplot(2,2,3)
hold on
%plot(finalTime,x(:,6),'*')
plot(finalTime, tangentialVelocity(:,3),'+')
ylabel('Computed Velocity (m/s)')
yyaxis right
plot(finalTime,groundTruthVelocity(:,3),'o')
hold off
title({'$v_z$ as a function of time'},'Interpreter','latex')
ylabel('Velocity (m/s)')
xlabel('Time (h)')
legend('Computed Tangential Velocity','Ground Truth Velocity')
% subplot(2,2,4)
% hold on
% plot(finalTime,vecnorm(x(:,4:6),2,1),'*')
% ylabel('Angular Velocity (rad/s)')
% yyaxis right
% plot(finalTime,vecnorm(groundTruthAng,2,1),'*')
% hold off
% title({'$\omega$ as a function of time'},'Interpreter','latex')
% ylabel('Ground Truth Angular Velocity (rad/s)')
% xlabel('Time (h)')

figure
hold on
plot(finalTime,x(:,4),'*')
plot(finalTime,x(:,5),'*')
plot(finalTime,x(:,6),'*')
plot(finalTime,groundTruthAng(:,1:3),'o')
hold off
title('Fitted versus Ground Truth Angular Velocities')
xlabel('Time (hours)')
ylabel('Angular Velocity m/s')
legend({'$\hat{\omega}_1$','$\hat{\omega}_2$','$\hat{\omega}_3$','$\omega_1$','$\omega_2$','$\omega_3$'},'Interpreter','latex')



% % we also need to propagate the quaternions since we need an estimate for
% % the querying of the star catalog. We're going to grab this from the
% % assignment we hade about the propagation of rotations.
% for i = 1:1%length(imds.Files)
%     [stars] = readStarCatalog(raDeg(1), decDeg(1), 0.85, 14.0);
%     [starStates, status] = getStarStates([betterQuaternion(i,:)' ,betterQuaternion(i,:)'], stars, [1024; 1024], 13*10^(-6), 0.85);
%     fprintf('Number of measured Stars is %d \n',size(dataVault(i).stars(:,1:2),1))
%     fprintf('Number of reference Stars is %d \n',size(starStates,2))
%     indexPointingStar(i).idx = star_pattern_matching_proj(starStates(1:2,:)',dataVault(i).stars(:,1:2));
%
%     trans = fitgeotrans(dataVault(1).stars(2:end,1:2),starStates(1:2,indexPointingStar(i).idx)','affine');
%     sss = trans.T(2, 1);
%     scc = trans.T(1, 1);
%     scale_recovered = sqrt(sss*sss+scc*scc);
%     translation_X_recovered = trans.T(3, 1);
%     translation_Y_recovered = trans.T(3, 2);
%     theta_recovered = atan2(sss, scc) * 180 / pi;
%
%     deltas = dataVault(1).stars(2:end,1:2) - starStates(1:2,indexPointingStar(i).idx)';
%     figure
%     plot(deltas(:,1),deltas(:,2),'*')


% %     figure
% %     hold on
% %     plot([dataVault(1).stars(:,1)], [dataVault(1).stars(:,2)],'*')
% %     plot(starStates(1,indexPointingStar(i).idx),starStates(2,indexPointingStar(i).idx),'*')
% %     hold off
% %
% %     figure
% %     hold on
% %     plot([dataVault(1).stars(:,1)], [dataVault(1).stars(:,2)],'*')
% %     plot(starStates(1,:),starStates(2,:),'*')
% %     hold off
% end












