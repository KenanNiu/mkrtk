clear 

initTestSuite;

hax = axes;
hold on
hp = patch();
hl = line;
hd = plot(0,0);

hax2 = axes;
hl2 = line;

% Invalid constructor calls
assertExceptionThrown( @()SelectiveBrush(),                 'SelectiveBrush:SelectiveBrush:IncorrectNargin' )
assertExceptionThrown( @()SelectiveBrush([hax,5,hl,hd,hp]), 'SelectiveBrush:SelectiveBrush:BadInputs' )
assertExceptionThrown( @()SelectiveBrush([hax,hl,hd,hp]),   'SelectiveBrush:SelectiveBrush:BadInputs' )
assertExceptionThrown( @()SelectiveBrush([hl,hl2]),         'SelectiveBrush:SelectiveBrush:ParentDiffers' )

% Valid constructor calls:
h = [hl,hd];
assert( isa( SelectiveBrush(hl) ,'SelectiveBrush') )
br = SelectiveBrush(hl);
br = SelectiveBrush(h);
assertEqual(br.axis, hax)
assertEqual(br.targets, h)


% ... There's more to be checked, but it didn't seem necessary at the time