% load original contact schedule
clc;
clear;
gait = 'spin_Laikago';
ctactSeq = readmatrix(sprintf('contact_%s.txt', gait)); % contact status of each foot at every time step
% k = 1;
% ctact_freq = [];
% l = 4;
% while k <= size(ctactSeq,1)
%     % if reach contact at k, look ahead
%     if ctactSeq(k, l)==1
%         j = 1;
%         while k+j <= size(ctactSeq,1)
%             if ~ctactSeq(k+j, l)
%                 ctact_freq = [ctact_freq, j];
%                 break;            
%             end
%             j = j + 1;
%         end
%         k = k + j + 1;
%     else
%         k = k + 1;
%     end
% end
% histogram(ctact_freq);

% ctactSeq(1:4,:) = ones(4,4);
% ctactSeq(5:6,:) = zeros(2,4);
% ctactSeq(11:13,:) = zeros(3,4);
% ctactSeq(14:19,:) = repmat([0,1,1,0],6,1);
% ctactSeq(20,:) = zeros(1,4);
% ctactSeq(21:26,:) = repmat([1,0,0,1],6,1);
% ctactSeq(27,:) = zeros(1,4);
% ctactSeq([29,30],:) = [0,1,1,0;0,1,1,0];
% ctactSeq(36, :) = [0,0,0,0];
% ctactSeq(37:42,:) = repmat([1,0,0,1], 6, 1);
% first_HR_ctact = 0;
% parttern = [repmat([0,1,1,0],6,1);
%             zeros(2,4);
%             repmat([1,0,0,1],6,1)];
% k = 42;
% while k <= size(ctactSeq,1)-13
%     if ctactSeq(k,3) == 0
%         first_HR_ctact = 1;
%         ctactSeq(k,:) = zeros(1,4);
%     end
%     if (first_HR_ctact) && (ctactSeq(k,3)==1)
%         first_HR_ctact = 0;
%         ctactSeq(k:k+13, :) = parttern;
%         k = k + 14;
%     else
%         k = k + 1;
%     end    
% end

% if gait == "pace_Laikago"
%     ctactSeq(:,1) = ctactSeq(:,3);
%     ctactSeq(:,2) = ctactSeq(:,4);
% end
if gait =="left_turn_A1"
    ctactSeq = filter_out_short_swing(ctactSeq, 6);
    ctactSeq(13:15, 4) = 1;
end
if gait == "spin_Laikago"
    ctactSeq = filter_out_short_swing(ctactSeq, 3);
end
len = 100;
X = categorical({'FR','FL','HR','HL'});
X = reordercats(X,{'FR','FL','HR','HL'});
Y = repmat(ones(1,len), 4, 1);
figure
b = barh(X, Y, 'stacked','FaceColor','flat');
for i=1:len
    for l = 1:4
        b(i).CData(l,:) = ones(1,3)-ctactSeq(i,l)*ones(1,3);
        drawnow
    end    
end
