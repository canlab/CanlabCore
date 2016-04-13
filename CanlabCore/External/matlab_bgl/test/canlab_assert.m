function canlab_assert(condition,varargin)

if condition, return;
else
  error(varargin{:});
end
