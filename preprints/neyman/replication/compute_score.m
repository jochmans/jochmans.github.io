%%% Computes the orthogonalized scores, for q=1 to 6, for computation of
%%% score differences in the Monte Carlo in the supplement
% Uses codes "team_CES_orderq_Teams_1_2_12" 

function obj=compute_score(par,Mat_meanoutput,Mat_production,Mat_alpha_prelim,qq)

Score=0;

for ii=1:size(Mat_meanoutput,1)
    if qq==1
        Score_i=team_CES_order1_Teams_1_2_12(Mat_meanoutput(ii,:),Mat_production(ii,:),...
            [par(1) par(2)],[par(4) par(3)],exp(Mat_alpha_prelim(ii,:)));
    elseif qq==2
        Score_i=team_CES_order2_Teams_1_2_12(Mat_meanoutput(ii,:),Mat_production(ii,:),...
            [par(1) par(2)],[par(4) par(3)],exp(Mat_alpha_prelim(ii,:)));
    elseif qq==3
        Score_i=team_CES_order3_Teams_1_2_12(Mat_meanoutput(ii,:),Mat_production(ii,:),...
            [par(1) par(2)],[par(4) par(3)],exp(Mat_alpha_prelim(ii,:)));
    elseif qq==4
        Score_i=team_CES_order4_Teams_1_2_12(Mat_meanoutput(ii,:),Mat_production(ii,:),...
            [par(1) par(2)],[par(4) par(3)],exp(Mat_alpha_prelim(ii,:)));
    elseif qq==5
        Score_i=team_CES_order5_Teams_1_2_12(Mat_meanoutput(ii,:),Mat_production(ii,:),...
            [par(1) par(2)],[par(4) par(3)],exp(Mat_alpha_prelim(ii,:)));
    elseif qq==6
        Score_i=team_CES_order6_Teams_1_2_12(Mat_meanoutput(ii,:),Mat_production(ii,:),...
            [par(1) par(2)],[par(4) par(3)],exp(Mat_alpha_prelim(ii,:)));
    end
    Score=Score+Score_i;
end

obj=Score;

end
