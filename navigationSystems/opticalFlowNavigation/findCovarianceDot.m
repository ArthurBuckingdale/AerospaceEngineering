function dot = findCovarianceDot(P_minus,F,G,W)
    Pdot = F*P_minus + P_minus*F' + G*W*G';
    dot = reshape(Pdot,[],1);
end