function seq = rsa_model_sequence(base_formula)
% rsa_model_sequence  Builder for nested LME model lists (sugar over rsa_compare_models).
%
% Builds a sequence of Wilkinson formulas by adding one term at a time.
% Drives rsa_compare_models without writing every formula by hand.
%
% Usage
% -----
%   seq = rsa_model_sequence('Y ~ 1 + (1|Sub)')
%   seq = seq.add_term('SameCondition');
%   seq = seq.add_term('SameBodysite');
%   seq = seq.add_term('SameCondition:SameBodysite');
%   seq = seq.add_term('SameSession');
%   seq = seq.upgrade_random('(SameCondition | Sub)');   % swap (1|Sub) -> (SameCondition|Sub)
%
%   formulas = seq.formulas;
%   [T, best] = dat.rsa_compare_models(formulas, ...);
%
% Inputs
% ------
%   base_formula  char/string Wilkinson formula to start from. Must contain
%                 a random-effects clause '(... | <var>)'.
%
% Returns
% -------
%   seq  struct with handles:
%     .add_term(term)         add a fixed-effect term to the current formula
%     .upgrade_random(re)     replace the random-effects clause
%     .formulas               cellstr of all formulas built so far
%     .current                most recent formula
%     .reset()                drop all but the base formula

if nargin < 1 || isempty(base_formula)
    error('rsa_model_sequence:noBase', 'Pass a base Wilkinson formula.');
end

% Internal state, captured by closures
state = struct();
state.base       = char(base_formula);
state.formulas   = {state.base};

seq = make_seq();

    function s = make_seq()
        s = struct();
        s.formulas        = state.formulas;
        s.current         = state.formulas{end};
        s.add_term        = @add_term_impl;
        s.upgrade_random  = @upgrade_random_impl;
        s.reset           = @reset_impl;
    end

    function s2 = add_term_impl(term)
        cur = state.formulas{end};
        new = inject_fixed_term(cur, char(term));
        state.formulas{end+1} = new;
        s2 = make_seq();
    end

    function s2 = upgrade_random_impl(new_re)
        cur = state.formulas{end};
        new = replace_random(cur, char(new_re));
        state.formulas{end+1} = new;
        s2 = make_seq();
    end

    function s2 = reset_impl()
        state.formulas = {state.base};
        s2 = make_seq();
    end

end


function out = inject_fixed_term(formula, term)
% Split off the random-effects clause, append the term to the fixed RHS,
% reattach the random clause.
[fixed_rhs, rand_clauses] = split_formula(formula);
% Add term if not already present
if ~ismember(strtrim(term), strsplit(strtrim(fixed_rhs), {'+', ' '}))
    if strcmp(strtrim(fixed_rhs), '1')
        new_fixed = term;
    else
        new_fixed = [strtrim(fixed_rhs) ' + ' term];
    end
else
    new_fixed = fixed_rhs;
end
out = ['Y ~ ' new_fixed];
if ~isempty(rand_clauses)
    out = [out ' + ' strjoin(rand_clauses, ' + ')];
end
end


function out = replace_random(formula, new_re)
[fixed_rhs, ~] = split_formula(formula);
out = ['Y ~ ' strtrim(fixed_rhs) ' + ' new_re];
end


function [fixed_rhs, rand_clauses] = split_formula(formula)
% Split a Wilkinson formula into fixed-effect RHS string and a list of
% random-effects parenthetical clauses.
rhs = regexprep(formula, '^\s*[A-Za-z0-9_]+\s*~\s*', '');
% Find parenthesized random-effect clauses
rand_clauses = regexp(rhs, '\([^)]*\)', 'match');
% Strip them out of the fixed RHS
fixed_rhs = regexprep(rhs, '\([^)]*\)', '');
% Clean up dangling +s and whitespace
fixed_rhs = regexprep(fixed_rhs, '\s*\+\s*\+\s*', ' + ');
fixed_rhs = regexprep(fixed_rhs, '^\s*\+\s*', '');
fixed_rhs = regexprep(fixed_rhs, '\s*\+\s*$', '');
fixed_rhs = strtrim(fixed_rhs);
if isempty(fixed_rhs), fixed_rhs = '1'; end
end
