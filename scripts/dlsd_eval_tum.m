%data = csvread('data/results.csv', 1, 0);
%mean(100*data(data(:, 2) > 0, 9)./data(data(:, 2) > 0, 2))
%mean(100*data(data(:, 4) > 0, 10)./data(data(:, 4) > 0, 4))
%mean(100*data(data(:, 6) > 0, 11)./data(data(:, 6) > 0, 6))

values = [1; 5; 10; 15; 20; 25; 30; 35; 40; 45; 50];
N = length(values);

lsd = zeros(N, 1);
ag3line = zeros(N, 1);
edlines = zeros(N, 1);

for idx = 1:N
    min_line_percent = values(idx);
    file = sprintf('data/results_%02d.csv', min_line_percent);
    data = csvread(file, 1, 0);
    lsd(idx) = mean(100*data(data(:, 2) > 0, 9)./data(data(:, 2) > 0, 2));
    ag3line(idx) = mean(100*data(data(:, 4) > 0, 10)./data(data(:, 4) > 0, 4));
    edlines(idx) = mean(100*data(data(:, 6) > 0, 11)./data(data(:, 6) > 0, 6));
    
    if values(idx) == 25
        data_25 = data;
    end
end

afigure(1);
hold on

plot(values, lsd)
plot(values, ag3line)
plot(values, edlines)

xlabel('Min. line length threshold [%]')
ylabel('Recall [%]')
legend('LSD', 'AG3lines', 'EDLines', 'Location', 'southeast')

disp("LSD:")
disp(sum(data_25(:, 2)))

disp("AG3line:")
disp(sum(data_25(:, 4)))

disp("EDLines:")
disp(sum(data_25(:, 6)))

disp("D-LSD:")
disp(sum(data_25(:, 8)))
