%% Part 0
close all;
clc;
clear;

%% Part 1
% Define our char array which contains alphabet and additional characters
chars = ['a':'z', ' ', '.', ',', '!', '"', ';'];
% chars_len = length(chars); % its 32

% Creating a mapset with 2 rows and columns : 
% Row 1 is the characters, Row 2 is the binary coded value of each character
Mapset = cell(2, 32);

% Set the values into Mapset
for i = 1:32
    % Store the character in the first row of Mapset
    Mapset{1, i} = chars(i);
    
    % Add the binary coded value for each character (i-1 in binary)
    Mapset{2, i} = dec2bin(i - 1, 5);
end

%% Part 3
message = 'signal';
coded_signal1 = coding_amp(message, 1, Mapset, 1);
coded_signal2 = coding_amp(message, 2, Mapset, 1);
coded_signal3 = coding_amp(message, 3, Mapset, 1);

%% Part 4 test
decoded_message1=decoding_amp(coded_signal1, 1, Mapset, 1);
decoded_message2=decoding_amp(coded_signal2, 2, Mapset, 1);
decoded_message3=decoding_amp(coded_signal3, 3, Mapset, 1);

disp(decoded_message1);
disp(decoded_message2);
disp(decoded_message3);

%% part 5
noise = randn(1,3000);

figure;
histogram(noise, 30, 'Normalization', 'pdf');
hold on;
x = linspace(min(noise), max(noise), 100);
plot(x, normpdf(x, 0, 1), 'r', 'LineWidth', 2);
title('Histogram and Gaussian Fit');
xlabel('Amplitude');
ylabel('Probability Density');
legend('Data Histogram', 'Gaussian Fit');


%% part 6
 
message = 'signal';
for i=1:3
    coded_signal = coding_amp(message, i, Mapset, 1);
    coded_signal = coded_signal + 0.01*randn(1, length(coded_signal));
    x_axis=zeros(1,length(coded_signal));
    for j=1:length(coded_signal)
        x_axis(j)=j/100;
    end
    figure;
    plot(x_axis,coded_signal)
    title('coded signal with noise');
    decoded_message=decoding_amp(coded_signal, i, Mapset, 1);
    
    disp([decoded_message, ' for BitRate = ', num2str(i)]);
end

%% part 7 and 8

message = 'signal';
for BitRate=1:3
    sum_var = 0;
    max_var = 0;
    for i=1:1000
        std = 0.01;
        decoded_message = 'signal';
        while message == decoded_message
            coded_signal = coding_amp(message, BitRate, Mapset, 0);
            coded_signal = coded_signal + std*randn(1, length(coded_signal));
            decoded_message=decoding_amp(coded_signal, BitRate, Mapset, 0);
            if message == decoded_message
                std = std + 0.01;
            end
        end
        sum_var = sum_var + (std^2);
        if max_var < (std-0.01)^2
            max_var =   (std-0.01)^2;
        end
    end
    x_axis=zeros(1,length(coded_signal));
    for j=1:length(coded_signal)
        x_axis(j)=j/100;
    end
    figure;
    plot(x_axis,coded_signal)
    title(['coded signal with noise with BitRate = ', num2str(BitRate) ]);
    disp([decoded_message, ' for BitRate = ', num2str(BitRate), ' mean noise variance = ', num2str(sum_var/1000)]);
    disp(['Max variance for BitRate = ', num2str(BitRate), ' decodes correctly is ', num2str(max_var)]);
end

%% part 11
message = 'signal';
for BitRate=1:3
    sum_var = 0;
    max_var = 0;
    for i=1:1000
        std = 0.01;
        decoded_message = 'signal';
        while message == decoded_message
            coded_signal = coding_amp(message, BitRate, Mapset, 0);
            coded_signal = coded_signal + std*randn(1, length(coded_signal));
            decoded_message=decoding_amp2(coded_signal, BitRate, Mapset, 0);
            if message == decoded_message
                std = std + 0.01;
            end
        end
        sum_var = sum_var + (std^2);
        if max_var < (std-0.01)^2
            max_var =   (std-0.01)^2;
        end
    end
    x_axis=zeros(1,length(coded_signal));
    for j=1:length(coded_signal)
        x_axis(j)=j/100;
    end
    figure;
    plot(x_axis,coded_signal)
    title(['coded signal with noise with BitRate = ', num2str(BitRate) ]);
    disp([decoded_message, ' for BitRate = ', num2str(BitRate), ' mean noise variance = ', num2str(sum_var/1000)]);
    disp(['Max variance for BitRate = ', num2str(BitRate), ' decodes correctly is ', num2str(max_var)]); 
end

%% part 2
function Coded_signal = coding_amp(message, BitRate, mapset, show_plt)

    fs = 100;  % Sampling frequency
    % Convert message to binary sequence using mapset
    message_bin = '';
    for k = 1:length(message)
        char_index = find(strcmp(mapset(1,:), message(k)));  % Find the character index in mapset
        if isempty(char_index)
            error(['Character ', message(k), ' not found in mapset.']);
        end
        char_bin = mapset{2, char_index};  % Get binary representation of each character
        message_bin = [message_bin, char_bin];  % Append to the binary message sequence
    end

    Coded_signal = [];
    message_bin_len = length(message_bin);

    t = 0.01:(1/fs):1;
    a = [];
    a = [a, zeros(1, length(t))];
    for i = 1:((2^BitRate)-1)
        a = [a, i/((2^BitRate)-1)*sin(2*pi*t)];
    end
    for i = 1:BitRate:message_bin_len
        code_number = bin2dec(message_bin(i:i+BitRate-1));
        Coded_signal = [Coded_signal, a((code_number)*100 + 1 : (code_number+1)*100)];
    end
    
    if (show_plt == 1)
        x_axis = zeros(1, message_bin_len/BitRate*100);
        for i = 1:message_bin_len/BitRate*100
            x_axis(i) = i/100;
        end
        figure;
        plot(x_axis, Coded_signal);
        grid on; 
        title('Coded Signal Generated from Message', 'FontSize', 14, 'FontWeight', 'bold');  
        xlabel('Time (s)', 'FontSize', 12);
        ylabel('Amplitude', 'FontSize', 12); 
        legend('Coded Signal', 'Location', 'best');  
    end
end

%% part 4
function decoded_message=decoding_amp(coded_signal, BitRate, mapset, show_plt)
    fs = 100;
    t = 0.01 : 0.01 : 1;
    time_of_message = (length(coded_signal)) / fs;

    correlation = zeros(1, time_of_message);

    Y = 2*sin(2*pi*t);
    for i = 1 : time_of_message
            X = coded_signal(fs*i-(fs-1):fs*i);
            correlation(i) = sum(X.*Y)/fs;
    end
    if show_plt == 1
        t = (0:length(correlation) - 1);
        figure;
        plot(t, correlation);
        title(['Convolution of Signal with One Period of 2*sin(2*pi*t) (bitrate = ' num2str(BitRate) ')']);
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
    % Find tresholds
    a = [];
    a = [a,0];
    for i=1:((2^BitRate)-1)
        a = [a, i/((2^BitRate)-1)];
    end
    treshold = [];
    for i=1:((2^BitRate)-1)
        treshold = [treshold, (a(i)+a(i+1))/2];
    end
    
    binarydecoded = [];
    for i=1:time_of_message
        if correlation(i) < treshold(1)
            binarydecoded = [binarydecoded dec2bin(0, BitRate)];
        end
        if correlation(i) > treshold((2^BitRate)-1)
            binarydecoded = [binarydecoded dec2bin(((2^BitRate)-1), BitRate)];
        end
    
        for j = 1 : ((2^BitRate)-1)-1
            if treshold(j) < correlation(i) && treshold(j+1) > correlation(i)
                binarydecoded = [binarydecoded dec2bin(j, BitRate)];
            end  
        end        
    end
    
    decoded_message = [];
    lenbin = length(binarydecoded);
    for i = 1 : lenbin/5 
        ch = binarydecoded(5*i-4 : 5*i);
        number = bin2dec(ch); 
        decoded_message = [decoded_message mapset{1, number+1}];
    end
end

%% part 11
function decoded_message=decoding_amp2(coded_signal, BitRate, mapset, show_plt)
    fs = 100;
    t = 0.01 : 0.01 : 1;
    time_of_message = (length(coded_signal)) / fs;

    correlation = zeros(1, time_of_message);

    Y = 10*sin(2*pi*t);
    for i = 1 : time_of_message
            X = coded_signal(fs*i-(fs-1):fs*i);
            correlation(i) = sum(X.*Y)/fs;
    end
    if show_plt == 1
        t = (0:length(correlation) - 1);
        figure;
        plot(t, correlation);
        title('Convolution of Signal with One Period of 2*sin(2*pi*t)');
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
    % Find tresholds
    a = [];
    a = [a,0];
    for i=1:((2^BitRate)-1)
        a = [a, 5*i/((2^BitRate)-1)];
    end
    treshold = [];
    for i=1:((2^BitRate)-1)
        treshold = [treshold, (a(i)+a(i+1))/2];
    end
    
    binarydecoded = [];
    for i=1:time_of_message
            if correlation(i) < treshold(1)
                binarydecoded = [binarydecoded dec2bin(0, BitRate)];
            end
        if correlation(i) > treshold((2^BitRate)-1)
            binarydecoded = [binarydecoded dec2bin(((2^BitRate)-1), BitRate)];
        end
    
        for j = 1 : ((2^BitRate)-1)-1
            if treshold(j) < correlation(i) && treshold(j+1) > correlation(i)
                binarydecoded = [binarydecoded dec2bin(j, BitRate)];
            end  
        end        
    end
    
    decoded_message = [];
    lenbin = length(binarydecoded);
    for i = 1 : lenbin/5 
        ch = binarydecoded(5*i-4 : 5*i);
        number = bin2dec(ch); 
        decoded_message = [decoded_message mapset{1, number+1}];
    end
end