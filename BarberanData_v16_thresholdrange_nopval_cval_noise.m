function [ n_nodes, n_edges, modular, cc, avgdeg, apl, diams ] = BarberanData_v16_thresholdrange_nopval_cval_noise( threshold_val, noise )
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Load data.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = load( 'soil454_otu_table_notext_noOTUname.txt' );
    data = data'; %tranpose the data
    [n, m] = size( data ); %n = number of locations, m = number of OTUs.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Add a column to index original locations.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loc_indices = 1:1:n;
    data = [ (1000.*loc_indices)' data ]; %multiply indices by 1000 to be sure 
    %they are the last column when we sort by column sums.


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Focusing on just OTUs with more than 5 sequences TOTAL:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    colsums = sum( data ); %total number of sequences in each OTU (total in each column)
    data5 = [ data ; colsums ]; %add totals as the last row so we can sort on that row.
    data5 = sortrows( data5', n+1 )'; %sort OTUs by the colsums value (number of sequences, ascending), in position 152 = n+1
    otucutoff = find( data5(n+1, :)>5, 1, 'first' ); %find the first OTU that has 6 or more sequences.
    chopdata = data5( 1:n, otucutoff:m+1 ); %get rid of colsums (last row), go from cutoff to m+1.


    %Move location codes from the last column back to the first column, and divide by 1000.
    [newn, newm] = size( chopdata );
    chopdata = [ (0.001.*chopdata( :, newm )), chopdata( :, 1:newm-1 ) ];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot of distribution of sequences per location, across all of the
    %locations.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % figure;
    % set(gca,'FontSize',14)
    % plot(sort(sum(chopdata(:,2:end),2)),'o')
    % xlabel('Sampling locations (sorted)')
    % ylabel('Sum of sequences at that location (minus rare OTUs)')
    % title('Distribution of number of sequences, total, per location')
    % 
    numreps = 10;

    n_nodes = zeros(1,numreps);
    n_edges = zeros(1,numreps);

    modular = zeros(1,numreps);
    cc = zeros(1,numreps);
    avgdeg = zeros(1,numreps);
    apl = zeros(1,numreps);
    diams = zeros(1,numreps);
    %cval = zeros(1,numreps);
    
    for z = 1:numreps
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This code is used to add noise and prevent tied ranks.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (noise == 1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %This code is used to add noise and prevent tied ranks.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            delta = min( diff( sort( unique( chopdata( :,2:end ) ) ) ) );
            %Random uniform values on the interval [a, b].
            a = -delta/1000;
            b = -a;
            epsilon = a + ( b-a ).*rand( newn,(newm-1) );
            noisy_data = [ chopdata( :,1 ), (chopdata( :,2:end ) + epsilon ) ];
        else
            noisy_data = [ chopdata( :,1 ), chopdata( :,2:end )];
        end;
    


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Here, we generate the "orderings" matrix.
        %We multiply data by -1 to sort locations by descending number of sequences.
        %But we leave location indices as is.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %noisy_data = [ noisy_data( :,1 ), -noisy_data( :,2:end ) ];
        %Sort the location column, based on each OTU column ("rankings").
        %First column of orderings is the location indices (which we need in order to compare OTUs next).
        order_indices = 1:1:newn;
        orderings = [ order_indices' zeros( newn, newm-1 ) ]; %initialize ordering array.
        
        if (noise == 1)
            for k = 2:newm
                %sort rows on that column, descending.
                %two columns: indices and kth column of neg. data.  Then, sort on second row, which is the kth column.
                temp1 = sortrows( [ noisy_data( :,1 ), noisy_data( :,k ) ] , -2); 
                temp2 = sortrows( [ temp1( :,1 ), order_indices'], 1 ); %put those ranks into the original order based on locs (1,2,3...151)
                orderings( :,k ) = temp2( :,2 ); %Now, write those ranks to the orderings file.
            end;
        else
            for k = 2:newm
                %two columns: indices and kth column of neg. data.  Then, sort on second row, which is the kth column.
                temp1 = sortrows( [ noisy_data( :,1 ), noisy_data( :,k ) ] , -2); 
                temp2 = tiedrank(temp1(:,2)); %get the tied ranks.
                temp3 = sortrows([temp1(:,1) temp2],1); %put those ranks into the original order based on locs (1,2,3...151)
                orderings(:,k) = temp3(:,2); %Now, write those ranks to the orderings file.
            end;
        end;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Here, we get the Spearman correlation coefficients and p-values.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [ rho, ~ ] = corr( orderings( :,2:end ),'Type','Spearman');

        for i=1:size(rho,1)
            rho(i,i)=0;
        end;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Extract the upper triangle of rho.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Prevent rho's = 0 from being removed by adding dummy value:
%         dummy = 0.000001;
%         rho = rho + dummy;
%         pval = pval + dummy;
%         %Next step: all non-upper-tri values set to zero.
%         rho_upper = triu(rho,1);
%         pval_upper = triu(pval,1);
%         %reshape into a vector.
%         rho_rs = reshape(rho_upper,1,(newm-1)^2);
%         pval_rs = reshape(pval_upper,1,(newm-1)^2);
%         %Then, remove all TRUE zeros (from lower tri):
%         rho_rs( :, all(~rho_rs,1)) = []; 
%         pval_rs( :, all(~pval_rs,1)) = [];
%         %finally, subtract the dummy value.
%         rho_rs = rho_rs - dummy;
%         pval_rs = pval_rs - dummy;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Extract the co-occurence network at rho > 0.6
        %We define the network as all co-occurrences with rho > 0.6 and pval < 0.01
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %network = ( abs( rho_upper ) > 0.6 ) & ( pval <= 0.01 );
        network = ( abs( rho ) > threshold_val );
        
        n_nodes(z) = sum( (sum( network ) > 0));
        n_edges(z) = sum( sum( network )) /2;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %AVERAGE DEGREE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        avgdeg(z) = (2*n_edges(z))/n_nodes(z);
        %cval(z) = (2*n_edges(z))/1577;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %DIAMETER
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dist = graphallshortestpaths(sparse(network));
        diams(z) = max(dist(~isinf(dist)));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %AVERAGE PATH LENGTH
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dist1 = dist(~isinf(dist));
        dist2 = dist1(dist1 ~=0);
        apl(z) = mean(dist2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Modularity. Edge list extracted and passed to R.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        network_upper_tri = triu(network,1);
        spar_n = sparse(network_upper_tri);
        nw_indices = find(spar_n); %this is a vector.
        edge_list = [mod(nw_indices,n_nodes(z)), ceil(nw_indices./n_nodes(z)) ];
        if (length(edge_list)) ~= 0
            for i = 1:length(edge_list(:,1))
                if edge_list(i,1)==0
                    edge_list(i,1)=n_nodes(z);
                end;
            end;
            %Next, I will pass edge-list to R to calculate the modularity.
            dlmwrite('edgelist.txt',edge_list,'delimiter',' ')
            unix('R CMD BATCH GetModularityGreedy.R');
            modularity = load('ModularityScore.txt');
            modular(z) = modularity;
        else modular(z) = NaN;
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Clustering coefficient.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        b = triu(network^2,1);
        c2 = sum(b); % number of connected triples
        c1 = sum(b.*network); % 3x number of triangles
        cc(z) = full(sum(c1)/sum(c2)); % clustering coefficient

    end;
    

    %Plot number of nodes
    %figure();
    %diff = max(n_nodes)-min(n_nodes);
    %y = hist(n_nodes, diff);
    %why = y./sum(y);
    %bar(min(n_nodes):max(n_nodes),why)
    %set(gca,'FontSize',14)
    %xlabel('Number of nodes')
    %ylabel('Frequency')
    %title('Distribution of number of nodes, 2000 reps, rho = 0.6, true data')

    %Plot number of edges
    %figure();
    %diff = max(n_edges)-min(n_edges);
    %x = hist(n_nodes,diff);
    %exe = x./sum(x);
    %bar(min(n_edges):max(n_edges),exe)
    %set(gca,'FontSize',14)
    %xlabel('Number of edges')
    %ylabel('Frequency')
    %title('Distribution of number of edges, 2000 reps, rho = 0.6, true data')



    %Plot Diameters 
    % figure();
    % diff = max(diams)-min(diams);
    % hist(diams,diff);
    % set(gca,'FontSize',14)
    % xlabel('Diameter')
    % ylabel('Frequency over 2000 trials')
    % hold on;
    % xlim([0 20])
    % title('Distribution of diameters, 2000 reps, rho=0.6, true data')
    % %line([18 18], [0.1 500],'Color',[1 0 0])

    %plot avg degrees
    % figure();
    % zyx = pdfsmooth(avgdeg');
    % plot(zyx(:,1),zyx(:,2),'LineWidth',2)
    % hold on;
    % set(gca,'FontSize',14)
    % xlabel('average degree')
    % ylabel('frequency over 2000 trials')
    % xlim([0 5])
    % title('Distribution of average degrees, 2000 reps, rho=0.6, true data')
    % line([4.59 4.59], [0 0.007],'Color',[1 0 0])

    %plot avg path length
    % wvu = pdfsmooth(apl');
    % plot(wvu(:,1),wvu(:,2), 'LineWidth',2)
    % xlim([1 6])
    % set(gca,'FontSize',14)
    % xlabel('Average Path Length')
    % ylabel('Frequency')
    % title('Distribution of avg. path lengths, 2000 reps, rho=0.6, true data')
    %line([5.53 5.53],[0 0.005],'Color',[1 0 0])


    %plot cc
    % tsr = pdfsmooth(cc');
    % plot(tsr(:,1),tsr(:,2),'LineWidth',2)
    % set(gca,'FontSize',14)
    % xlabel('Clustering Coefficient')
    % ylabel('Frequency')
    % title('Distribution of clustering coeff., 2000 reps, rho=0.6, true data')
    % xlim([0 1])
    %line([0.33 0.33],[0 0.003],'Color',[1 0 0])

    %plot modular
    % qpo = pdfsmooth(modular');
    % plot(qpo(:,1),qpo(:,2),'LineWidth',2)
    % set(gca,'FontSize',14)
    % xlabel('Modularity')
    % ylabel('Frequency')
    % xlim([0 1])
    % title('Distribution of modularity scores, 2000 reps, rho=0.6, true data')
    %line([0.77 0.77],[0 0.003],'Color',[1 0 0])

end
