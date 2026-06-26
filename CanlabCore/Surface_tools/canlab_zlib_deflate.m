function comp = canlab_zlib_deflate(raw)
% Deflate (compress) bytes to a raw zlib stream natively via the JVM.
%
% :Usage:
% ::
%     comp = canlab_zlib_deflate(raw)
%
% Produces a raw zlib stream (header 0x78 0x9C) suitable for GIFTI/CIFTI
% "GZipBase64Binary" payloads after base64 encoding. Uses java.util.zip
% (ships with MATLAB) -- no external toolbox.
%
% :Inputs:
%   **raw:** uint8 vector of bytes to compress.
%
% :Outputs:
%   **comp:** uint8 column vector, the raw zlib-compressed stream.
%
% :See also: canlab_zlib_inflate, canlab_write_gifti, canlab_write_cifti

raw = uint8(raw(:));
baos = java.io.ByteArrayOutputStream();
dos  = java.util.zip.DeflaterOutputStream(baos);
dos.write(typecast(raw, 'int8'), 0, numel(raw));
dos.finish();
dos.close();
comp = typecast(int8(baos.toByteArray()), 'uint8');
comp = comp(:);
end
