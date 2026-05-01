function obj = get_sample_fmri_data()
%GET_SAMPLE_FMRI_DATA Load the standard emotionreg sample for tests.
%
%   obj = get_sample_fmri_data() returns the 30-image Wager 2008 emotion
%   regulation contrast set as an fmri_data object. Used by tests to keep
%   load semantics in one place; if load_image_set's behavior shifts, only
%   this helper needs updating.

obj = load_image_set('emotionreg', 'noverbose');
end
