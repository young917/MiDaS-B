for dataname={"threads-ask-ubuntu", "coauth-MAG-Geolgoy-full", "coauth-MAG-History-full"};
for dataname={"coauth-MAG-History-full"};
    for portion_str={"0.1", "0.3", "0.5", "0.7", "0.9"};
        [fileID, msg] = fopen('input/answer/' + dataname{1} +  '/icd_row_' + portion_str{1} + '.txt','r');
        if (fileID < 0);
            disp(msg);
            continue;
        end;
        fileID2 = fopen('input/answer/' + dataname{1} +  '/icd_col_' + portion_str{1} + '.txt','r');
        fileID3 = fopen('input/answer/' + dataname{1} +  '/dim_' + portion_str{1} + '.txt','r');
        formatSpec = '%f';
        R = fscanf(fileID, formatSpec);                       
        C = fscanf(fileID2, formatSpec);
        Dim = fscanf(fileID3, formatSpec);
        S = sparse(R,C,1,Dim(1), Dim(2)); 
        sv = svds(S, 300);
        sz = size(sv,1);
        dlmwrite('output/answer/' + dataname{1} +  '/sv_full_' + portion_str{1} + '.txt',sv(:),'newline','pc');
    end;
end;