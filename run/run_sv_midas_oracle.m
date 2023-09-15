for repeatname={"1", "2", "3"};
    for dataname={"threads-ask-ubuntu", "coauth-MAG-Geolgoy-full", "coauth-MAG-History-full"};
        for portion_str={"0.1", "0.3", "0.5", "0.7", "0.9"};
            for alphaname={"0.0000", "0.2500", "0.5000", "1.0000", "2.0000"};
                [fileID, msg] = fopen('input/es/add_global_deg_min_' + alphaname{1} + '/' + dataname{1} +  '/' + repeatname{1} + '/icd_row_' + portion_str{1} + '.txt','r');
                if (fileID < 0);
                    disp(msg);
                    continue;
                end;
                fileID2 = fopen('input/es/add_global_deg_min_' + alphaname{1} + '/' + dataname{1} +  '/' + repeatname{1} + '/icd_col_' + portion_str{1} + '.txt','r');
                fileID3 = fopen('input/es/add_global_deg_min_' + alphaname{1} + '/' + repeatname{1} + '/dim_' + portion_str{1} + '.txt','r');
                formatSpec = '%f';
                R = fscanf(fileID, formatSpec);                       
                C = fscanf(fileID2, formatSpec);
                Dim = fscanf(fileID3, formatSpec);
                S = sparse(R,C,1,Dim(1), Dim(2)); 
                sv = svds(S, 300);
                sz = size(sv,1);
                dlmwrite('output/es/add_global_deg_min_' + alphaname{1} + '/' + dataname{1} +  '/' + repeatname{1} + '/sv_full_' + portion_str{1} + '.txt',sv(:),'newline','pc');
            end;
        end;
    end;
end;