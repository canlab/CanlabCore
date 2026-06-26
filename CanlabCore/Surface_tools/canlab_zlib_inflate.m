function out = canlab_zlib_inflate(raw)
% Inflate (decompress) a raw zlib byte stream natively via the JVM.
%
% :Usage:
% ::
%     out = canlab_zlib_inflate(raw)
%
% GIFTI "GZipBase64Binary" and similar payloads are raw zlib streams (header
% 0x78 0x9C), NOT gzip. This decompresses a uint8 vector entirely on the JVM
% side, which avoids a long-standing MATLAB bug where an in-place
% java.util.zip.Inflater.inflate(buffer) call returns garbage because MATLAB
% passes the Java method a copy of the buffer. Uses java.util.zip plus
% org.apache.commons.io.IOUtils (both ship with MATLAB) -- no external toolbox.
%
% :Inputs:
%   **raw:** uint8 vector, a raw zlib-compressed stream.
%
% :Outputs:
%   **out:** uint8 column vector of the decompressed bytes.
%
% :See also: canlab_zlib_deflate, canlab_read_gifti, canlab_read_cifti

bais = java.io.ByteArrayInputStream(typecast(uint8(raw(:)), 'int8'));
iis  = java.util.zip.InflaterInputStream(bais);
baos = java.io.ByteArrayOutputStream();
org.apache.commons.io.IOUtils.copy(iis, baos);
iis.close();
out = typecast(int8(baos.toByteArray()), 'uint8');
out = out(:);
end
