close all; clear;

numreps = 1000;
dim = 300;
% nodes_matrix = zeros(number_of_rand_inputs,numreps);
% edges_matrix = zeros(number_of_rand_inputs,numreps);
% mod_matrix = zeros(number_of_rand_inputs,numreps);
% cc_matrix = zeros(number_of_rand_inputs,numreps);
% diam_matrix = zeros(number_of_rand_inputs,numreps);
% degree_matrix = zeros(number_of_rand_inputs,numreps);
% path_matrix = zeros(number_of_rand_inputs,numreps);

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data = load('soil454_otu_table_notext_noOTUname.txt');
% data = data'; %tranpose the data
% [n, m] = size( data ); %n = number of locations, m = number of OTUs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add a column to index original locations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loc_indices = 1:1:n;
% data = [ (1000.*loc_indices)' data ]; %multiply indices by 1000 to be sure 
% %they are the last column when we sort by column sums.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Focusing on just OTUs with more than 5 sequences TOTAL:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colsums = sum( data ); %total number of sequences in each OTU (total in each column)
% data5 = [ data ; colsums ]; %add totals as the last row so we can sort on that row.
% data5 = sortrows( data5', n+1 )'; %sort OTUs by the colsums value (number of sequences, ascending), in position 152 = n+1
% otucutoff = find( data5(n+1, :)>5, 1, 'first' ); %find the first OTU that has 6 or more sequences.
% chopdata = data5( 1:n, otucutoff:m+1 ); %get rid of colsums (last row), go from cutoff to m+1.
% 
% 
% %Move location codes from the last column back to the first column, and divide by 1000.
% [newn, newm] = size( chopdata );
% chopdata1 = [ (0.001.*chopdata( :, newm )), chopdata( :, 1:newm-1 ) ];

%These steps are new for Method 2:
%chopdata2 = chopdata1(:,(2:end));
%chopdata_reshape = reshape(chopdata2, 1, newn*(newm-1));


%for eye = 1:number_of_rand_inputs

    %rand_data = randsample(chopdata_reshape, newn*(newm-1), false);
    %rand_data_rs = reshape(rand_data,newn,newm-1);
    %chopdata3 = [ loc_indices', rand_data_rs];

%Initialize vectors for the repetitions:
n_nodes = ones(1,numreps);
n_edges = zeros(1,numreps);

modular = zeros(1,numreps);
cc = zeros(1,numreps);
avgdeg = zeros(1,numreps);
avgdeg_includingzeros = zeros(1,numreps);
apl = zeros(1,numreps);
diams = zeros(1,numreps);

for z = 1:numreps

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This code is used to add noise and prevent tied ranks.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     delta = min( diff( sort( unique( chopdata3( :,2:end ) ) ) ) );
%     a = -delta/1000;
%     b = -a;
%     epsilon = a + ( b-a ).*rand( newn,(newm-1) );
%     noisy_data = [ chopdata3( :,1 ), (chopdata3( :,2:end ) + epsilon ) ];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Here, we generate the "orderings" matrix.
    %We multiply data by -1 to sort locations by descending number of sequences.
    %But we leave location indices as is.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     noisy_data = [ noisy_data( :,1 ), noisy_data( :,2:end ) ];
%     %Sort the location column, based on each OTU column ("rankings").
%     %First column of orderings is the location indices (which we need in order to compare OTUs next).
%     order_indices = 1:1:newn;
%     orderings = [ order_indices' zeros( newn, newm-1 ) ]; %initialize ordering array.

%     for k = 2:newm
%         %sort rows on that column, descending.
%         %two columns: indices and kth column of neg. data.  Then, sort on second row, which is the kth column.
%         temp1 = sortrows( [ noisy_data( :,1 ), noisy_data( :,k ) ] , -2); 
%         temp2 = sortrows( [ temp1( :,1 ), order_indices'], 1 ); %put those ranks into the original order based on locs (1,2,3...151)
%         orderings( :,k ) = temp2( :,2 ); %Now, write those ranks to the orderings file.
%     end;
% 
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Here, we get the Spearman correlation coefficients and p-values.
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      [ rho, pval ] = corr( orderings( :,2:end ),'Type','Spearman');



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract the upper triangle of rho.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %Prevent rho's = 0 from being removed by adding dummy value:
%     dummy = 0.000001;
%     rho = rho + dummy;
%     pval = pval + dummy;
%     %Next step: all non-upper-tri values set to zero.
%     rho_upper = triu(rho,1);
%     pval_upper = triu(pval,1);
%     %reshape into a vector.
%     rho_rs = reshape(rho_upper,1,(newm-1)^2);
%     pval_rs = reshape(pval_upper,1,(newm-1)^2);
%     %Then, remove all TRUE zeros (from lower tri):
%     rho_rs( :, all(~rho_rs,1)) = []; 
%     pval_rs( :, all(~pval_rs,1)) = [];
%     %finally, subtract the dummy value.
%     rho_rs = rho_rs - dummy;
%     pval_rs = pval_rs - dummy;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Extract the co-occurence network at rho > 0.6
    %We define the network as all co-occurrences with rho > 0.6 and pval < 0.01
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %network = ( abs( rho_upper ) > 0.6 ) & ( pval <= 0.01 );
    %network = ( abs( rho ) > 0.36); %& ( pval <= 0.01 );

    %%%%%%%%
    %I believe we just put in a new Erdos-Renyi graph, right here
    %%%%%%%%


    network = rand(dim,dim);
    network_upper = triu(network,1);
    network = network_upper>0.9504; %this threshold 1746.07 edges for 300 nodes, to get the right avg deg.
    for i = 1:dim
        for j = i+1:dim
            if network(i,j) == 1
                network(j,i) =1;
            end;
        end;
    end
    
%     
%     b = +network;
%     issymmetric(b);
  
    
    n_nodes(z) = sum( (sum( network ) > 0));
    n_edges(z) = sum( sum( network ))/2;


%         figure;
%         hist(rho_rs,100);
%         xlim([-1 1]);
%         set(gca,'FontSize',14)
%         title('Distribution of rho in pairwise Spearman analysis of 454 data')
%         xlabel('Rho')
%         ylabel('Frequency')



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %AVERAGE DEGREE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    avgdeg(z) = n_edges(z)/n_nodes(z);
    %avgdeg(z) = (2*n_edges(z))/n_nodes(z);
    %avgdeg_includingzeros(z) = (2*n_edges(z))/dim;

   %{
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
%}
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
%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Clustering coefficient.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b = triu(network^2,1);
    c2 = sum(b); % number of connected triples
    c1 = sum(b.*network); % 3x number of triangles
    cc(z) = full(sum(c1)/sum(c2)); % clustering coefficient




    %Plot distribution of rho and pval; 
    % %rho_rs = reshape(rho,1,(newm-1)^2);
    % pval_rs = reshape(pval,1,(newm-1)^2);
    % 
    % figure;
    % hist(rho_rs,100);
    % xlim([-1 1]);
    % set(gca,'FontSize',14)
    % title('Distribution of rho in pairwise Spearman analysis of 454 data')
    % xlabel('Rho')
    % ylabel('Frequency')

    %Alternatively, to use PDFSmooth rather than hist:
    %g = pdfsmooth(rho_rs');
    %plot(g(:,1),g(:,2))

    % 
    % figure;
    % hist(pval_rs,100);
    % xlim([0 1]);
    % set(gca,'FontSize',14)
    % title('Distribution of p-values in pairwise Spearman of 454 data')
    % xlabel('p-value')
    % ylabel('Frequency')
    % 
    % 
    % figure;
    % hist(pval_rs,1000);
    % xlim([0 0.1]);


    % Remove zero rows
%     network( all(~network,2), : ) = [];
% 
%     % Remove zero columns
%     network( :, all(~network,1) ) = [];

    %networkcolsums = sum(network);
    %network1 = [network; networkcolsums];
    %[rows cols] = size(network1);
    %sortnetw = sortrows(network1',rows)';
    %networkCO=find(sortnetw(rows,:)>0,1,'first');
    %nw_minus_singletons = sortnetw(networkCO:rows-1,:);


    %sumagain = sum(nw_minus_singletons);
    %network2 = [nw_minus_singletons; sumagain];
    %[rows cols] = size(network2);
    %sortnetwagain = sortrows(network2',rows)';
    %networkcutoff = find(sortnetwagain(rows,:)>0,1,'first');


%}
end;
    
% % nodes_matrix(eye,:) = n_nodes(:);
% % edges_matrix(eye,:) = n_edges(:);
% % mod_matrix(eye,:) = modular(:);
% % cc_matrix(eye,:) = cc(:);
% % degree_matrix(eye,:) = avgdeg(:);
% % path_matrix(eye,:) = apl(:);
% % diam_matrix(eye,:) = diams(:);
%     
%     %Plot Diameters 
% figure();
% hist(diams);
% xlabel('Diameter')
% ylabel('Frequency over 2000 trials')
% hold on;
% xlim([0 20])
% line([18 18], [0.1 500],'Color',[1 0 0])
% 
%     %plot avg degrees
% figure();
% zyx = pdfsmooth(avgdeg');
% plot(zyx(:,1),zyx(:,2),'LineWidth',2)
% hold on;
% xlabel('average degree')
% ylabel('frequency over 2000 trials')
% xlim([0 5])
% line([4.59 4.59], [0 0.007],'Color',[1 0 0])
% 
%     %plot avg path length
% figure();
% wvu = pdfsmooth(apl');
% plot(wvu(:,1),wvu(:,2), 'LineWidth',2)
% hold on;
% xlim([1 6])
% line([5.53 5.53],[0 0.005],'Color',[1 0 0])
% xlabel('Average Path Length')
% ylabel('Frequency')
% 
% %plot cc
% figure();
% tsr = pdfsmooth(cc');
% plot(tsr(:,1),tsr(:,2),'LineWidth',2)
% line([0.33 0.33],[0 0.003],'Color',[1 0 0])
% xlabel('Clustering coefficient')
% ylabel('Frequency')
    % close;
%end;

%csvwrite('ErdosRenyiMethod_apl.csv',apl)
%csvwrite('ErdosRenyiMethod_cc.csv',cc)
%csvwrite('ErdosRenyiMethod_diams.csv',diams)
%csvwrite('ErdosRenyiMethod_deg.csv',avgdeg)
%csvwrite('ErdosRenyiMethod_mod.csv',modular)

toc
