function convmodel = plot_models_and_hrf(model, modelname, dt, scalehrf)

%% Model 1

if ~exist('dt', 'var')
    dt = 0.0606;
end 

% HRF
hrf = spm_hrf(dt);

if exist('scalehrf', 'var')
    convmodel = conv(model, hrf *scalehrf);
else
    convmodel = conv(model, hrf);
    % Scale convmodel to plot it in the same order of magnitude of model
    k = max(model)/max(convmodel);% scaling factor
    convmodel = convmodel * k;
end

% convmodel_adapted = convmodel(1:size(model, 2));
% model = [model zeros(1, 528)]




t = 0:dt:90;
t2 = 0:dt:122;
figure
hold on
p1 = plot(t, model, 'blue')
p2 = plot(t2, convmodel, 'red', 'LineWidth', 2)
% p1 = plot(model, 'blue')
% p2 = plot(convmodel, 'red')

xlabel('Time (s)')
ylabel('Modelled activity')


if exist('modelname', 'var')
    title(modelname)
end

grid on

xline(0, ':', 'cue', 'LabelVerticalAlignment', 'bottom');
xline(2, ':', 'trial onset', 'LabelVerticalAlignment', 'bottom');
xline(54, ':', 'trial end', 'LabelVerticalAlignment', 'bottom');
xline(56, ':', 'resting state', 'LabelVerticalAlignment', 'bottom');
xline(90, ':', 'end block', 'LabelVerticalAlignment', 'bottom');

limmax = max([model, convmodel])+ 2;
limmin = min([model, convmodel])- 2;
axis([0 122 limmin limmax]);

legend([p1 p2], 'model', 'conv(model, hrf)')
hold off

end