DataX(Gencount,:,irun) = BestParamX;
DataY(Gencount,:,irun) = BestParamY;
DataF(irun) = CostF;
if CnstrFlag == 0
    DataF(irun) = 1e10;
end

numColumn = 2;
if done == true
    [~, minindex] = min(DataF);
    fprintf('Plotting Data of Run = %d\n', minindex);
    figure;
    set(gcf,'position',[0,0,200 * numColumn,200 * ceil(size(DataX,2)/numColumn)])
    for j=1:size(DataX,2)
        subplot(ceil(size(DataX,2)/numColumn),numColumn,j);
        plot(1:Gencount, DataX(:,j,minindex));
        if j > size(DataX,2) - numColumn
            xlabel('Generation');
        end
        ylabel(sprintf('x%d',j));
        %axis([0 Gencount ParamMinX(j) ParamMaxX(j)]);
        grid on;
    end
    print(gcf,sprintf('results/%s_x_%s.png', CostDef, method),'-dpng','-r300')

    figure;
    set(gcf,'position',[200 * numColumn,0,200 * numColumn,200 * ceil(size(DataY,2)/numColumn)])
    for j=1:size(DataY,2)
        subplot(ceil(size(DataY,2)/numColumn),numColumn,j);
        plot(1:Gencount, DataY(:,j,minindex));
        if j > size(DataY,2) - numColumn
            xlabel('Generation');
        end
        ylabel(sprintf('y%d',j));
        %axis([0 Gencount ParamMinX(j) ParamMaxX(j)]);
        grid on;
    end
    print(gcf,sprintf('results/%s_y_%s.png', CostDef, method),'-dpng','-r300')
    pause(1);
end