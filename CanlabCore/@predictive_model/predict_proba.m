function P = predict_proba(obj, X)
% predict_proba  Calibrated probability of the positive class.
%
% Requires a prior calibrate() call (which stashes the Platt /
% isotonic calibrator in obj.fitted_values.calibrator).
%
% :Usage:
% ::
%     pm = calibrate(pm, X, Y);
%     P = predict_proba(pm, X_new);   % [n x 1] in [0,1]
%
% :Inputs:
%
%   **obj:**
%        a calibrated @predictive_model (run calibrate() first).
%
%   **X:**
%        [n x p] predictor matrix in the model's feature space.
%
% :Outputs:
%
%   **P:**
%        [n x 1] calibrated probability of the positive class, in [0, 1].
%
% :Examples:
% ::
%     dat = load_image_set('DPSP_hotwarm');
%     X = dat.dat'; Y = dat.Y;
%     pm = predictive_model('algorithm','svm','task','classification');
%     pm = calibrate(pm, X, Y, 'method', 'platt');
%     P  = predict_proba(pm, X);       % P(class = +1)
%
% :See also:
%   calibrate, predict, fit

    if ~isfield(obj.fitted_values, 'calibrator') || isempty(obj.fitted_values.calibrator)
        error('predictive_model:predict_proba:NotCalibrated', ...
            'No calibrator on object. Call calibrate(pm, X, Y) first.');
    end
    cal = obj.fitted_values.calibrator;

    [~, scores] = predict(obj, X);
    s = scores(:, end);

    switch cal.method
        case 'platt'
            P = 1 ./ (1 + exp(-(cal.A * s + cal.B)));
        case 'isotonic'
            % Interpolate the PAVA solution at s; clamp to the
            % extremes of the calibration set.
            P = interp1(cal.s_knots, cal.p_knots, s, 'previous', NaN);
            below = s < cal.s_knots(1);
            above = s > cal.s_knots(end);
            P(below) = cal.p_knots(1);
            P(above) = cal.p_knots(end);
        otherwise
            error('predictive_model:predict_proba:UnknownCalibrator', ...
                'Unknown calibrator method: %s', cal.method);
    end
    P = max(min(P, 1), 0);
end
