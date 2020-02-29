clear all;
close all;
%% read in the image
im = imread('test3.jpg');
%% color space transformation
% change color space to hsv
hsv = rgb2hsv(im); 
% change color space to lab
lab = rgb2lab(im);
[H, L, ~] = size(lab);
% statistic information in 2d hist
temp = lab(:,:,2);
a = temp(:);
temp = lab(:,:,3);
b = temp(:);
h = hsv(:,:,1);
%% build up groud truth for this photo
% need to be edit mannully per sample
gnd_true =false(H, L);
% gnd_true(382:386,819:823) = 1; %% label test1.jpg
% gnd_true(917:921,862:866) = 1;
% gnd_true(919:923,1555:1559) = 1;
% gnd_true(515:515+8,783:783+8) = 1; %% label test2.jpg
% gnd_true(1378:1378+8,3160:3160+8) = 1;
% gnd_true(1594:1594+8,1752:1752+8) = 1;
gnd_true(1219:1219+4,399:399+4) = 1; %% label test3.jpg
gnd_true(260:264,536:540) = 1;
gnd_true(728:728+4,1561:1565) = 1;
figure('Name','gnd_t')
imshow(gnd_true);
%% knn search
X = [a,b];
Y = X;
kn = 26; % just above 25
[Idx, D] = knnsearch(X,Y,'K',kn);
P = kn./(length(a).*(pi.*(D(:,kn)+10e-4)).^2); %% +samll number incase 0
P = reshape(P,H,L);
[counts, centers] = hist(D(:,kn), 100);
% how to decide the distance
sum_1 = 0;
cent = 0;
for i = fliplr(1: length(counts))
    if counts(i) < 30 % this is an emperical number, which means the upper 
        % bound of area of anomalous object. 100? 60? 30?
        % test1 30, test2 60
        sum_1 = sum_1 + counts(i);
    else
        cent = i;
        break;
    end
%     if counts(i) > 5000000
%         cent = i;
%         break
%     end
end
label = D(:,kn) > centers(cent); % decide by distance
label = reshape(label,H,L);
% show the result of knn
figure('Name','knn')
%% improvment for False Alarm
% CC = bwconncomp(label);
% numPixels = cellfun(@numel,CC.PixelIdxList);
% decision = numPixels > 1;
% for i = 1 : length(decision)
%     if (decision(i))
%         label(CC.PixelIdxList{i}) = 0;
%     end
% end
%% find the center of each object before MRF
% center = regionprops(label,'centroid');
% centroids = round(cat(1,center.Centroid)); % round to int
% % swap the coordinate
% temp_1 = centroids(:,1);
% centroids(:,1) = centroids(:,2);
% centroids(:,2) = temp_1;
% [num_ano,~] = size(centroids);
% %% calculate (inner/outer)GMM for each interest area
% discre1 = [];
% for i = 1: num_ano
%     ano = [];
%     nor = [];
%     try
%         for j  = centroids(i,1)-7:centroids(i,1)+7
%             for k = centroids(i,2)-7:centroids(i,2)+7
%                 if (label(j,k) == 1)
%                     ano = cat(2,ano,h(j,k));
%                 else
%                     nor = cat(2,nor,h(j,k));
%                 end
%             end
%         end
%     catch
%         try
%             centroids(i,:)=[];
%             label(centroids(i,1)-7:centroids(i,1),centroids(i,2)-7:centroids(i,2)) = 0;
%         catch
%         end
%         continue
%     end
%     discre1 = cat(1, discre1,[myGMM(0:0.0001:1,ano,0.0006),myGMM(0:0.0001:1,nor,0.0006)]);
% end
imshow(label);
%% find the threshold for probability
% [row, column] = find(label);
% em = [];
% for i = 1: length(row)
%     em = [em,P(row(i),column(i))];
% end
% thresh_p = max(em);
% label_p = P <= thresh_p; % decide threshold for p
% %% the part of foreground
% P1 = P;
% mu = [0 0];
% sigma = [6 0; 0 6];
% m_mask = ones(7);
% m_mask(4,4) = 0; % create a mask to count
% for n = 1:5
% for i = 1:length(row) % detect on each anomalious
%     row(i)
%     column(i)
%     if row(i)-3<0 || row(i)+3>H || column(i)-3 < 0 || column(i)+3 > L
%         continue; % incase outrange
%     end
%     
%     temp_list = [];
%     for j = row(i)-1: row(i)+1 % original is 3, change to 1
%         for k = column(i)-1: column(i)+1 % span a 7*7 window to detect
%             if j == row(i) && k == column(i) % don't compare itself
%                 continue;
%             end
%             if label_p(j,k) == 1
%                 temp_1 = lab(j,k,2:3)-lab(row(i),column(i),2:3);
%                 temp_1 = reshape(temp_1,1,2);
%                 temp_list = [temp_list, ...
%                   sqrt((2*pi)^2*36)*mvnpdf(temp_1,mu,sigma)];
%               % sqrt((2*pi)^2*36)*
%             end
%         end
%     end
%     if isempty(temp_list)
%         continue; % which means no anomalious in neighbourhood
%     else 
%         P1(row(i),column(i)) = P1(row(i),column(i))/mean(temp_list);
%     end
% end
% label_p1 = P1 <= thresh_p;
% end
% figure('Name','foreground')
% 
% imshow(label_p1)
%%             
n_mask = ones(3);
n_mask(2,2) = 0; % create a mask to perform Ising potential
label_1 = cell(1,4);
label_1{1} = label;
gama = 0.5*10^(-4); % 1*10-4
T = 0.75;

for i =1:2
    anom_normal_cnt = -8+2.*conv2(double(label_1{i}),n_mask,'same'); % Qf-Qb 
    % count the number of anomalous 
    if i ==1 %% for parameter gama adjustment
        P_t = gama.*exp((1/T).*anom_normal_cnt);
    end
    label_1{i+1} = P<gama.*exp((1/T).*anom_normal_cnt); % improb is the things unchange
end
figure('Name','knn with MRF')
imshow(label_1{3})
%% find the center of each object after MRF
final_out = label_1{3};
center = regionprops(final_out,'centroid');
centroids = round(cat(1,center.Centroid)); % round to int
% swap the coordinate
temp_1 = centroids(:,1);
centroids(:,1) = centroids(:,2);
centroids(:,2) = temp_1;
[num_ano,~] = size(centroids);
%% calculate (inner/outer)GMM for each interest area
% window size is 7+7+1  15*15
discre2 = [];
for i = 1: num_ano
    ano = [];
    nor = [];
    try
        for j  = centroids(i,1)-7:centroids(i,1)+7
            for k = centroids(i,2)-7:centroids(i,2)+7
                try
                if (label(j,k) == 1)
                    ano = cat(2,ano,h(j,k));
                else
                    nor = cat(2,nor,h(j,k));
                end
                catch
                    continue
                end
            end
        end
    catch
    end
    discre2 = cat(1,discre2,[myGMM(0:0.0001:1,ano,0.004),myGMM(0:0.0001:1,nor,0.004)]); 
    % test1  0.004     test2 0.05
end
%% eliminate FP
[t,q] = size(discre2);
dist = zeros(t,1);
for i = 1 : t
    dist(i) = KLDiv(discre2(i,1:q/2),discre2(i,q/2+1:end));
end

[idx,C] = kmeans(dist,2);
if C(1) > C(2)
    TP = 1;
else
    TP = 2;
end
CC = bwconncomp(final_out);
for i = 1: length(idx)
    if (idx(i) ~= TP)
        final_out(CC.PixelIdxList{i}) = 0;
    end
end
figure('Name','knn with MRF after GMM and KLdiv')
imshow(final_out)





%%

% scatter(a, b, '*');
% title('Scatter of a & b');
% xlabel('a');
% ylabel('b');

% hist(D(:,1)-D(:,26))
