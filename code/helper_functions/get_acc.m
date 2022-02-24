% Compute classification accuracy for 3D motion decoding experiment

function accuracy = get_acc(predictions, p)
    test_labels = round(p.stimval(p.test_trials)*8/180);
    new_est = round(predictions/(180/8));
    new_est(new_est==0)=8; % make sure labels are [1:8]
    %accuracy = ones(1, length(test_labels));
    %accuracy(test_labels ~= new_est) = 0;
    %mean(accuracy);
    accuracy = mean(test_labels == new_est);
end